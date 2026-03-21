#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TattleTail: A Pyocin Prediction Tool (V1.0) 

Features:
- BLAST database self-check
- Mixed FFN/FASTA input
- Prodigal prediction or direct translation
- BLASTP search, functional label mapping, window-based clustering, 
- Four-step candidate evaluation
- Batch summary reports (plain text + colored version)
- Cluster hits annotations
 
Per-sample output files:
  - query_proteins.faa
  - blast_results.txt
  - tailocin_report.txt
  - tailocin_report.json
  - run.log
  - cluster_annotation.tsv

Batch summary files generated in root directory (plain + colored versions)
"""


import os
import re
import sys
import argparse
import subprocess
import json
from shutil import which
from datetime import datetime
from collections import defaultdict

from termcolor import colored
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------------------------- ANSI tools ---------------------------

ANSI_GREEN = "\033[32m"
ANSI_RED   = "\033[31m"
ANSI_RESET = "\033[0m"

def colorize_pass_fail_block(text: str) -> str:
    text = re.sub(r"\bPASS\b", f"{ANSI_GREEN}PASS{ANSI_RESET}", text)
    text = re.sub(r"\bFAIL\b", f"{ANSI_RED}FAIL{ANSI_RESET}", text)
    text = re.sub(r"(verdict:\s*)(TAILOCIN CANDIDATE ✅)", r"\1" + ANSI_GREEN + r"\2" + ANSI_RESET, text, flags=re.IGNORECASE)
    text = re.sub(r"(verdict:\s*)(not a candidate ❌)",   r"\1" + ANSI_RED   + r"\2" + ANSI_RESET, text, flags=re.IGNORECASE)
    text = re.sub(r"(Conclusion:\s*)(1 .*identified\.)", r"\1" + ANSI_GREEN + r"\2" + ANSI_RESET, text)
    text = re.sub(r"(Conclusion:\s*)(0 .*identified\.)", r"\1" + ANSI_RED   + r"\2" + ANSI_RESET, text)
    return text

def open_log_append(path):
    return open(path, "a") if path else open(os.devnull, "w")

def split_csv(s):
    return [x.strip() for x in s.split(",") if x.strip()]

# --------------------------- NEW: BLAST DB 自检 ---------------------------

def ensure_blast_db(db_fasta="database_tailocin.fasta", db_name="database_tailocin"):
    """
    Check if the BLAST protein database is ready and consistent.
    Exits with error message and makeblastdb command if not.
    Args:
        db_fasta: Path to the FASTA file used to build the database
                  (also used for mapping subject ID → functional label)
        db_name:  Prefix name of the BLAST database (used in blastp -db)
    """
    # 1) Check if the FASTA exists
    if not os.path.exists(db_fasta):
        if os.path.basename(db_fasta) == "database_tailocin.fasta":
            print("database_tailocin.fasta not exists", file=sys.stderr)
        else:
            print(f"[FATAL] DB FASTA not exists: {db_fasta}", file=sys.stderr)
        sys.exit(2)

    # 2) Check if the index file exists: (protein library:. bin/. phr/. psq)
    idx_files = [db_name + ".pin", db_name + ".phr", db_name + ".psq"]
    have_index = all(os.path.exists(p) for p in idx_files)
    if not have_index:
        if which("makeblastdb") is None:
            print("[FATAL] makeblastdb not exists", file=sys.stderr)
            print("        please install BLAST+: (conda install -c bioconda blast)", file=sys.stderr)
        cmd = f'makeblastdb -in "{db_fasta}" -dbtype prot -out "{db_name}" -parse_seqids'
        print("[FATAL] BLAST database index not found:", db_name, file=sys.stderr)
        print("        Please create an index first (copy the next command to run):", file=sys.stderr)
        print("        " + cmd, file=sys.stderr)
        sys.exit(2)

    # 3) Attempt to process:
    ok = True
    if which("blastdbcmd") is not None:
        try:
            r = subprocess.run(["blastdbcmd", "-db", db_name, "-info"],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if r.returncode != 0:
                ok = False
        except Exception:
            ok = False
    if not ok:
        cmd = f'makeblastdb -in "{db_fasta}" -dbtype prot -out "{db_name}" -parse_seqids'
        print("[FATAL] Database index found to exist but unreadable:", db_name, file=sys.stderr)
        print("        Suggest rebuilding the index:", file=sys.stderr)
        print("        " + cmd, file=sys.stderr)
        sys.exit(2)

    # 4) FASTA is newer than index file: WARNING
    try:
        fasta_mtime = os.path.getmtime(db_fasta)
        idx_mtime = min(os.path.getmtime(p) for p in idx_files)
        if fasta_mtime > idx_mtime:
            cmd = f'makeblastdb -in "{db_fasta}" -dbtype prot -out "{db_name}" -parse_seqids'
            print("[WARN] DB FASTA is newer than BLAST index. Suggest rebuilding the index.", file=sys.stderr)
            print("       " + cmd, file=sys.stderr)
    except Exception:
        pass

    # 5) Quick examination header format (at least two tokens, the second being a functional tag)
    bad = 0
    try:
        with open(db_fasta, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().split()
                    if len(parts) < 2:
                        bad += 1
        if bad > 0:
            print(f"[FATAL] In {db_fasta} there has {bad} header(s) didn't include second functional tag", file=sys.stderr)
            print("        Must meet the format: '>ID LABEL [optional]'", file=sys.stderr)
            print("        Example: '>BAR70105.1 integrase_1 [Pseudomonas aeruginosa]'", file=sys.stderr)
            sys.exit(2)
    except Exception:
        pass

# --------------------------- Dependency and Input Types ---------------------------

def check_dependencies(require_prodigal=True, require_blastp=True):
    ok = True
    if require_prodigal and which("prodigal") is None:
        print("prodigal not exists", file=sys.stderr); ok = False
    if require_blastp and which("blastp") is None:
        print("blastp not exists", file=sys.stderr); ok = False
    if not ok:
        sys.exit(1)

def is_annotated_fasta(file_path):
    annotated_keywords = ["gene=", "cds", "product=", "locus_tag", "[location"]
    with open(file_path, 'r', errors="ignore") as f:
        data = f.read(1024)
    return any(k in data.lower() for k in annotated_keywords)

# --------------------------- FFN→FAA Translate ---------------------------

def translate_ffn_to_faa(ffn_file, faa_file, genetic_code=11, partial_policy="trim", log_path=None):
    def log(msg):
        if log_path:
            with open(log_path, "a") as lf:
                lf.write("[FFN_TRANSLATE] " + msg + "\n")

    records = []
    bad_chars = re.compile(r"[^ACGTN]")

    for record in SeqIO.parse(ffn_file, "fasta"):
        raw = str(record.seq).upper()
        clean = raw.replace("U", "T").replace("-", "")
        clean = re.sub(r"\s+", "", clean)
        if bad_chars.search(clean):
            new_clean = re.sub(bad_chars, "", clean)
            if new_clean != clean:
                log(f"{record.id}: removed non-ACGTN chars; {len(clean)}->{len(new_clean)} bp")
            clean = new_clean

        rem = len(clean) % 3
        if rem != 0:
            if partial_policy == "trim":
                log(f"{record.id}: length {len(clean)} not multiple of 3; trimming tail {rem} nt")
                clean = clean[:len(clean) - rem]
            elif partial_policy == "pad":
                pad = "N" * (3 - rem)
                log(f"{record.id}: length {len(clean)} not multiple of 3; padding '{pad}'")
                clean = clean + pad
            elif partial_policy == "skip":
                log(f"{record.id}: length {len(clean)} not multiple of 3; skipped")
                continue
            elif partial_policy == "error":
                raise ValueError(f"{record.id}: CDS length {len(clean)} not multiple of 3")
            else:
                log(f"{record.id}: unknown partial_policy={partial_policy}; trimming by default")
                clean = clean[:len(clean) - rem]

        prot = Seq(clean).translate(table=genetic_code, to_stop=True)
        protein_record = SeqRecord(prot, id=record.id, description=f"{record.description} | cds_len={len(clean)} policy={partial_policy}")
        records.append(protein_record)

    SeqIO.write(records, faa_file, "fasta")

# --------------------------- Prodigal ---------------------------

def prodigal_mode_fix(mode):
    mode = str(mode).lower()
    return "meta" if mode == "meta" else "single"

def run_prodigal(nucl_fasta, faa_out, prodigal_mode="meta", genetic_code=11, log_path=None):
    cmd = ["prodigal", "-i", nucl_fasta, "-a", faa_out, "-p", prodigal_mode_fix(prodigal_mode), "-q", "-g", str(genetic_code)]
    try:
        with open_log_append(log_path) as lf:
            lf.write(f"[CMD] {' '.join(cmd)}\n")
            subprocess.run(cmd, check=True, stdout=lf, stderr=lf)
        return True, ""
    except subprocess.CalledProcessError as e:
        msg = f"Prodigal failed: {e}"
        if log_path:
            with open_log_append(log_path) as lf:
                lf.write(f"[ERROR] {msg}\n")
        return False, msg

# --------------------------- BLASTP ---------------------------

def run_blastp(query_faa, db_name, out_file, evalue="1e-10", threads=None, log_path=None):
    cmd = [
        "blastp", "-query", query_faa,
        "-db", db_name, "-out", out_file,
        "-evalue", str(evalue),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    ]
    if threads and int(threads) > 0:
        cmd += ["-num_threads", str(int(threads))]
    try:
        with open_log_append(log_path) as lf:
            lf.write(f"[CMD] {' '.join(cmd)}\n")
            subprocess.run(cmd, check=True, stdout=lf, stderr=lf)
        return True, ""
    except subprocess.CalledProcessError as e:
        msg = f"BLASTP failed: {e}"
        if log_path:
            with open_log_append(log_path) as lf:
                lf.write(f"[ERROR] {msg}\n")
        return False, msg

# --------------------------- Functional Mapping ---------------------------

def parse_functional_mapping(database_fasta):
    id_to_function = {}
    id_to_label = {}
    with open(database_fasta, errors="ignore") as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split()
                prot_id = parts[0][1:]
                label = parts[1]  # e.g., integrase_1
                function = label.split('_')[0].lower()
                id_to_function[prot_id] = function
                id_to_label[prot_id] = label
    return id_to_function, id_to_label

def parse_all_hits(blast_file, id_to_function, id_to_label,
                   evalue_thr=1e-10, pident_thr=35.0, length_thr=50, bitscore_thr=50.0):
    hits = []
    found_functions = set()
    if not os.path.exists(blast_file):
        return hits, found_functions

    with open(blast_file, errors="ignore") as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            qseqid = cols[0]; sseqid = cols[1]
            try:
                pident = float(cols[2]); alen = int(cols[3])
                evalue = float(cols[10]); bits = float(cols[11])
            except ValueError:
                continue

            if evalue > evalue_thr:   continue
            if pident < pident_thr:   continue
            if alen   < length_thr:   continue
            if bits   < bitscore_thr: continue

            func  = id_to_function.get(sseqid, "unknown")
            label = id_to_label.get(sseqid, sseqid)
            if func == "unknown":
                continue

            hits.append({
                "q": qseqid, "s": sseqid,
                "func": func, "label": label,
                "pident": pident, "length": alen,
                "bitscore": bits, "evalue": evalue
            })
            found_functions.add(func)
    return hits, found_functions

# --------------------------- Parse Prodigal Coords ---------------------------

def parse_prodigal_coords(faa_file):
    qcoords = {}
    for rec in SeqIO.parse(faa_file, "fasta"):
        hdr = rec.description
        try:
            parts = hdr.split('#')
            if len(parts) >= 3:
                left  = parts[0].strip()
                start = int(parts[1].strip())
                end   = int(parts[2].strip())

                base = left.split()[0].lstrip('>')
                contig = base
                if "_" in base:
                    prefix, tail = base.rsplit("_", 1)
                    if tail.isdigit():
                        contig = prefix

                qcoords[rec.id] = (contig, min(start, end), max(start, end))
        except Exception:
            continue
    return qcoords

def parse_ffn_header_coords(faa_file, log_path=None):
    """
    Recognized formats include:
      >lcl|NC_004547.2_cds_... [location=complement(32..475)]
      >ref|NZ_CP014999.1_cds_... [location=<1179..>2171]
      >... [location=join(100..200,300..400)]
    Return: { qseqid: (contig, start, end) }
    """
    def log(msg):
        if log_path:
            with open(log_path, "a") as lf:
                lf.write("[FFN_COORDS] " + msg + "\n")

    qcoords = {}
    for rec in SeqIO.parse(faa_file, "fasta"):
        hdr = rec.description
        token = hdr.split()[0]  # 如 lcl|NC_004547.2_cds_WP_...

        # grab <CONTIG> from <prefix>|<CONTIG>_cds...  
        m = re.search(r'^[^|]+\|(.+?)_cds', token)
        if m:
            contig = m.group(1)
        else:
            m2 = re.search(r'^[^|]+\|([^\s]+)', token)
            contig = m2.group(1) if m2 else re.split(r'_cds', token)[0]

        # extract [location=...]
        m = re.search(r'\[location=([^\]]+)\]', hdr)
        if not m:
            continue
        loc = m.group(1)

        # trim complement()/join()/order()
        loc_clean = loc
        for kw in ['complement', 'join', 'order']:
            loc_clean = re.sub(rf'{kw}\(', '(', loc_clean)
        loc_clean = loc_clean.replace('(', '').replace(')', '')

        starts, ends = [], []
        for part in loc_clean.split(','):
            part = part.strip()
            if not part:
                continue
            # Allow </> prefix
            mm = re.match(r'[<>]?(\d+)\.\.[<>]?(\d+)', part)
            if mm:
                s = int(mm.group(1)); e = int(mm.group(2))
            else:
                mm2 = re.match(r'[<>]?(\d+)$', part)
                if not mm2:
                    continue
                s = e = int(mm2.group(1))
            starts.append(s); ends.append(e)

        if not starts:
            continue
        s = min(starts); e = max(ends)
        if s > e:
            s, e = e, s

        qcoords[rec.id] = (contig, s, e)

    if qcoords:
        log(f"Parsed {len(qcoords)} coordinate entries from FFN-style headers.")
    else:
        log("No FFN-style coordinates parsed.")
    return qcoords

# --------------------------- Clustering & Evaluation ---------------------------

def join_hits_with_coords(hits, qcoords):
    ghits = []
    for h in hits:
        q = h["q"]
        if q in qcoords:
            contig, s, e = qcoords[q]
            gh = dict(h)
            gh["contig"] = contig
            gh["start"]  = s
            gh["end"]    = e
            ghits.append(gh)
    return ghits

def cluster_by_window(ghits, window_bp=10000):
    by_contig = defaultdict(list)
    for gh in ghits:
        by_contig[gh["contig"]].append(gh)

    clusters = []
    for contig, arr in by_contig.items():
        arr.sort(key=lambda x: x["start"])
        cur = []; cur_s = None; cur_e = None
        for h in arr:
            s, e = h["start"], h["end"]
            if not cur:
                cur = [h]; cur_s, cur_e = s, e
            else:
                if s - cur_e <= window_bp:
                    cur.append(h); cur_e = max(cur_e, e)
                else:
                    clusters.append({"contig": contig, "start": cur_s, "end": cur_e, "members": cur})
                    cur = [h]; cur_s, cur_e = s, e
        if cur:
            clusters.append({"contig": contig, "start": cur_s, "end": cur_e, "members": cur})
    return clusters

def evaluate_clusters_strict(
    clusters,
    flank=("trpE","trpG"),
    blockers=("capsid","terminase","integrase"),
    regulators=("prtN","prtR"),
    toxins=("holin","endolysin")
):
    results = []
    set_flank = set(x.lower() for x in flank)
    set_block = set(x.lower() for x in blockers)
    set_reg   = set(x.lower() for x in regulators)
    set_tox   = set(x.lower() for x in toxins)

    for clu in clusters:
        funcs = set(m["func"].lower() for m in clu["members"])

        step1_ok = set_flank.issubset(funcs)
        step2_ok = (len(funcs & set_block) == 0)
        step3_ok = set_reg.issubset(funcs)
        step4_ok = set_tox.issubset(funcs)

        results.append({
            "contig": clu["contig"],
            "start":  clu["start"],
            "end":    clu["end"],
            "members": clu["members"],
            "funcs":   funcs,
            "step1_ok": step1_ok, "step1_reason": None if step1_ok else "Missing flanking pair (trpE & trpG)",
            "step2_ok": step2_ok, "step2_reason": None if step2_ok else "Blocker(s) present (capsid/terminase/integrase)",
            "step3_ok": step3_ok, "step3_reason": None if step3_ok else "Missing regulators (prtN & prtR)",
            "step4_ok": step4_ok, "step4_reason": None if step4_ok else "Missing lysis genes (holin & endolysin)",
            "is_candidate": (step1_ok and step2_ok and step3_ok and step4_ok)
        })
    return results

def trim_cluster_by_flanks(cluster_eval, flank_markers=("trpE","trpG")):
    """
    Redefine cluster boundaries as: after trpE end to before trpG start.
    Handles reversed order by sorting flanks by position.
    Remove flank genes from members.
    Only applied if cluster is_candidate == True and both flanks exist.
    """
    flank_set = set(x.lower() for x in flank_markers)

    for c in cluster_eval:
        if not c["is_candidate"]:
            continue

        # Find trpE & trpG members
        flank_hits = [m for m in c["members"] if m["func"].lower() in flank_set]

        if len(flank_hits) < 2:
            # log: missing flanks
            print(f"[TRIM WARN] Cluster {c['contig']}:{c['start']}-{c['end']} missing both flanks for trim; skipping.", file=sys.stderr)
            continue

        trpE_hit = next((m for m in flank_hits if m["func"].lower() == "trpe"), None)
        trpG_hit = next((m for m in flank_hits if m["func"].lower() == "trpg"), None)

        if not trpE_hit or not trpG_hit:
            print(f"[TRIM WARN] Cluster {c['contig']}:{c['start']}-{c['end']} missing trpE or trpG; skipping trim.", file=sys.stderr)
            continue

        # Log: flanks location
        print(f"[TRIM INFO] Cluster {c['contig']}:{c['start']}-{c['end']}", file=sys.stderr)
        print(f"  trpE: {trpE_hit['start']}-{trpE_hit['end']}", file=sys.stderr)
        print(f"  trpG: {trpG_hit['start']}-{trpG_hit['end']}", file=sys.stderr)

        # Sort by start, find the left and right flags (ignore tags)
        flanks_sorted = sorted([trpE_hit, trpG_hit], key=lambda x: x["start"])
        left_flank = flanks_sorted[0]
        right_flank = flanks_sorted[1]

        # Calculate the core boundary: from left end+1 to right start -1
        new_start = left_flank["end"] + 1
        new_end = right_flank["start"] - 1

        # Check negative span
        if new_start > new_end:
            print(f"[TRIM WARN] Negative core span: {new_start} > {new_end}; possible assembly issue or non-adjacent flanks. Skipping trim.", file=sys.stderr)
            continue  # 跳过 trim，保持原 cluster，但仍移除 flanks 成员

        # update cluster coords
        c["start"] = new_start
        c["end"] = new_end

        # remove all flanks from members, and filter out genes that exceed the new boundary
        new_members = []
        for m in c["members"]:
            if m["func"].lower() in flank_set:
                continue
            if m["start"] >= new_start and m["end"] <= new_end:
                new_members.append(m)

        c["members"] = new_members

        # Update functions 
        c["funcs"] = set(m["func"].lower() for m in new_members)

        # Output log
        print(f"  Computed core: {new_start}..{new_end} ({(new_end - new_start + 1)/1000:.2f} kb)", file=sys.stderr)
        print(f"  Retained members: {len(new_members)}", file=sys.stderr)

    return cluster_eval




def top_hits_by_function(hits):
    best = {}
    for h in hits:
        f = h["func"].lower()
        if (f not in best) or (h["bitscore"] > best[f]["bitscore"]):
            best[f] = {
                "query": h["q"],
                "label": h["label"],
                "pident": h["pident"],
                "length": h["length"],
                "bitscore": h["bitscore"],
                "evalue": h["evalue"],
                "contig": h.get("contig"),
                "start":  h.get("start"),
                "end":    h.get("end")
            }
    return best

# --------------------------- Output ---------------------------
def write_report_and_json(
    report_txt, report_json, sample_name, params,
    found_functions_global, clusters_eval,
    flank_markers, regulatory_markers, toxin_markers, phage_blockers,
    block_scope, binary_output
):
    n_cand = sum(1 for c in clusters_eval if c["is_candidate"])
    if binary_output:
        conclusion_text = ("1 potential tailocin has been identified."
                           if n_cand > 0
                           else "0 potential tailocin-encoding biosynthetic gene clusters were identified.")
    else:
        conclusion_text = f"{n_cand} candidate cluster(s) identified."

    with open(report_txt, "w", encoding="utf-8") as rpt:
        rpt.write("Tailocin Finder Report (cluster-aware)\n")
        rpt.write("=====================================\n\n")
        rpt.write(f"Sample: {sample_name}\n")
        rpt.write(f"Params: {json.dumps(params, ensure_ascii=False)}\n\n")

        rpt.write("Note: Candidate clusters have been trimmed to exclude flanking trpE and trpG genes.\n")
        rpt.write("      The coordinates and gene list below represent the core pyocin region only.\n\n")

        rpt.write("Global summary (across all hits):\n")
        def _line(label, keys):
            states = [f"{k}:{'Y' if k.lower() in found_functions_global else 'N'}" for k in keys]
            return f"  {label}: " + ", ".join(states) + "\n"
        rpt.write(_line("Flanks", flank_markers))
        rpt.write(_line("Regulators", regulatory_markers))
        rpt.write(_line("Toxins", toxin_markers))
        rpt.write(_line("Blockers", phage_blockers))
        rpt.write(f"  Block scope: {block_scope}\n\n")

        if clusters_eval:
            rpt.write("Cluster-level results:\n")
            for idx, c in enumerate(clusters_eval, 1):
                span_kb = (c["end"] - c["start"] + 1) / 1000.0
                rpt.write(f"Cluster #{idx}  contig={c['contig']}  core span (trimmed)={c['start']}..{c['end']}  ({span_kb:.2f} kb)\n")
                funcs_sorted = sorted(c["funcs"])
                rpt.write(f"  functions (core pyocin genes): {', '.join(funcs_sorted)}\n")
                rpt.write("  Steps:\n")
                rpt.write(f"    1) flanks (trpE & trpG) detected outside core: {'PASS' if c['step1_ok'] else 'FAIL'}"
                          + ("" if c['step1_ok'] else f" — {c['step1_reason']}") + "\n")
                rpt.write(f"    2) blockers absent: {'PASS' if c['step2_ok'] else 'FAIL'}"
                          + ("" if c['step2_ok'] else f" — {c['step2_reason']}") + "\n")
                rpt.write(f"    3) regulators (prtN & prtR): {'PASS' if c['step3_ok'] else 'FAIL'}"
                          + ("" if c['step3_ok'] else f" — {c['step3_reason']}") + "\n")
                rpt.write(f"    4) lysis (holin & endolysin): {'PASS' if c['step4_ok'] else 'FAIL'}"
                          + ("" if c['step4_ok'] else f" — {c['step4_reason']}") + "\n")
                verdict = "TAILOCIN CANDIDATE ✅" if c["is_candidate"] else "not a candidate ❌"
                rpt.write(f"  verdict: {verdict}\n")

                if c["is_candidate"]:
                    rpt.write("\n  Core genes in this candidate cluster (sorted by genomic position):\n\n")
                    rpt.write("  function    label           locus_tag                   contig               start       end      %id   len   bitscore   evalue\n")
                    rpt.write("  " + "-" * 110 + "\n")

                    members_sorted = sorted(c["members"], key=lambda x: x["start"])
                    for m in members_sorted:
                        contig_disp = m["contig"][:20] + "..." if len(m["contig"]) > 20 else m["contig"]
                        locus_disp  = m["q"][:24] + "..." if len(m["q"]) > 24 else m["q"]
                        rpt.write(
                            f"  {m['func']:<11} {m['label']:<15} {locus_disp:<27} {contig_disp:<20} "
                            f"{m['start']:>9} {m['end']:>9} {m['pident']:>5.1f} {m['length']:>5} "
                            f"{m['bitscore']:>9.1f} {m['evalue']:>9.2e}\n"
                        )
                    rpt.write("\n")
                else:
                    rpt.write("\n")

        else:
            rpt.write("No clusters produced or no hits with coordinates (see run.log).\n\n")

        rpt.write(f"Conclusion: {conclusion_text}\n")

    # JSON section remain unchanged
    out = {
        "sample": sample_name,
        "params": params,
        "block_scope": block_scope,
        "binary_output": bool(binary_output),
        "binary_conclusion": 1 if (n_cand > 0) else 0,
        "candidates": int(n_cand),
        "global_found_functions": sorted(list(found_functions_global)),
        "clusters": []
    }
    for idx, c in enumerate(clusters_eval, 1):
        cluster_dict = {
            "cluster_id": idx,
            "contig": c["contig"],
            "start": c["start"],
            "end": c["end"],
            "functions": sorted(list(c["funcs"])),
            "is_candidate": c["is_candidate"],
            "steps": {
                "step1_flanks": {"ok": c["step1_ok"], "reason": c["step1_reason"]},
                "step2_blockers_absent": {"ok": c["step2_ok"], "reason": c["step2_reason"]},
                "step3_regulators": {"ok": c["step3_ok"], "reason": c["step3_reason"]},
                "step4_lysis": {"ok": c["step4_ok"], "reason": c["step4_reason"]}
            },
            "top_hits": top_hits_by_function(c["members"])
        }
        # Optional: If you want to record the member list in JSON as well
        if c["is_candidate"]:
            cluster_dict["core_members"] = [
                {
                    "function": m["func"],
                    "label": m["label"],
                    "locus_tag": m["q"],
                    "contig": m["contig"],
                    "start": m["start"],
                    "end": m["end"],
                    "pident": m["pident"],
                    "length": m["length"],
                    "bitscore": m["bitscore"],
                    "evalue": m["evalue"]
                }
                for m in sorted(c["members"], key=lambda x: x["start"])
            ]
        out["clusters"].append(cluster_dict)

    with open(report_json, "w", encoding="utf-8") as jf:
        json.dump(out, jf, indent=2, ensure_ascii=False)

def print_report_to_console(sample_name, params, found_functions_global, clusters_eval,
                            flank_markers, regulatory_markers, toxin_markers, phage_blockers,
                            block_scope, binary_output):
    n_cand = sum(1 for c in clusters_eval if c["is_candidate"])
    if binary_output:
        conclusion_text = ("1 potential tailocin has been identified."
                           if n_cand > 0
                           else "0 potential tailocin-encoding biosynthetic gene clusters were identified.")
    else:
        conclusion_text = f"{n_cand} candidate cluster(s) identified."

    print("Tailocin Finder Report (cluster-aware)")
    print("=====================================\n")
    print(f"Sample: {sample_name}")
    print(f"Params: {json.dumps(params, ensure_ascii=False)}\n")

    print(colored("Note: Candidate clusters trimmed — excluding trpE & trpG flanks", "yellow"))
    print("      Shown below: core pyocin region only.\n")

    print("Global summary (across all hits):")
    def yn(key): return 'Y' if key.lower() in found_functions_global else 'N'
    print("  Flanks:     " + ", ".join([f"{k}:{yn(k)}" for k in flank_markers]))
    print("  Regulators: " + ", ".join([f"{k}:{yn(k)}" for k in regulatory_markers]))
    print("  Toxins:     " + ", ".join([f"{k}:{yn(k)}" for k in toxin_markers]))
    print("  Blockers:   " + ", ".join([f"{k}:{yn(k)}" for k in phage_blockers]))
    print(f"  Block scope: {block_scope}\n")

    if clusters_eval:
        print("Cluster-level results:")
        for idx, c in enumerate(clusters_eval, 1):
            span_kb = (c["end"] - c["start"] + 1) / 1000.0
            print(f"Cluster #{idx}  contig={c['contig']}  core span (trimmed)={c['start']}..{c['end']}  ({span_kb:.2f} kb)")
            funcs_sorted = sorted(c["funcs"])
            print(f"  functions (core): {', '.join(funcs_sorted)}")
            print("  Steps:")
            def pf(ok): return colored("PASS", "green") if ok else colored("FAIL", "red")
            print(f"    1) flanks detected outside: {pf(c['step1_ok'])}" + ("" if c['step1_ok'] else f" — {c['step1_reason']}"))
            print(f"    2) blockers absent:          {pf(c['step2_ok'])}" + ("" if c['step2_ok'] else f" — {c['step2_reason']}"))
            print(f"    3) regulators:               {pf(c['step3_ok'])}" + ("" if c['step3_ok'] else f" — {c['step3_reason']}"))
            print(f"    4) lysis genes:              {pf(c['step4_ok'])}" + ("" if c['step4_ok'] else f" — {c['step4_reason']}"))
            verdict = colored("TAILOCIN CANDIDATE ✅", "green") if c["is_candidate"] else colored("not a candidate ❌", "red")
            print(f"  verdict: {verdict}")

            if c["is_candidate"]:
                print("\n  Core genes (sorted by position):")
                print(colored("  function    label           locus_tag                   contig               start       end      %id   len   bitscore   evalue", "cyan"))
                print("  " + "-" * 110)

                members_sorted = sorted(c["members"], key=lambda x: x["start"])
                for m in members_sorted:
                    contig_disp = m["contig"][:20] + "..." if len(m["contig"]) > 20 else m["contig"]
                    locus_disp  = m["q"][:24] + "..." if len(m["q"]) > 24 else m["q"]
                    print(
                        f"  {m['func']:<11} {m['label']:<15} {locus_disp:<27} {contig_disp:<20} "
                        f"{m['start']:>9} {m['end']:>9} {m['pident']:>5.1f} {m['length']:>5} "
                        f"{m['bitscore']:>9.1f} {m['evalue']:>9.2e}"
                    )
                print("")
            else:
                print("")

    else:
        print("No clusters produced or no hits with coordinates (see run.log).\n")

    print(f"Conclusion: {conclusion_text}")


# --------------------------- Batch compile ---------------------------
def aggregate_batch_text(results_root, items, timestamp):
    n = len(items)
    base = f"allresultsof{n}samples_{timestamp}"
    plain_path = os.path.join(results_root, base + ".txt")
    colored_path = os.path.join(results_root, base + "_colored.txt")
    os.makedirs(results_root, exist_ok=True)

    # ====================== PLAIN TEXT version ======================
    with open(plain_path, "w", encoding="utf-8") as out:
        out.write(f"# Tailocin Finder Batch Report — {n} samples — {timestamp}\n")
        out.write("# Each section below is the full content of per-sample tailocin_report.txt.\n\n")
        
        out.write("## SUMMARY\n")
        # New header: Add candidate_stans and cluster_stezes_kb
        out.write("No.\tsample\tcandidates\tbinary\tcandidate_spans\tcluster_sizes_kb\treport_dir\n")
        
        for i, it in enumerate(items, 1):
            spans = it.get("candidate_spans", "")
            sizes = it.get("cluster_sizes_kb", "")
            out.write(f"{i}\t{it['sample']}\t{it['n_cand']}\t{it['binary']}\t{spans}\t{sizes}\t{it['dir']}\n")
        
        out.write("\n## DETAILS\n\n")
        
        for i, it in enumerate(items, 1):
            out.write(f"### [{i}] Sample: {it['sample']}  ({it['dir']})\n\n")
            try:
                with open(it["report_txt"], "r", encoding="utf-8") as rf:
                    content = rf.read()
                out.write(content.rstrip() + "\n\n")
            except Exception as e:
                out.write(f"[ERROR] Cannot read {it['report_txt']}: {e}\n\n")

    # ====================== COLORED version ======================
    with open(colored_path, "w", encoding="utf-8") as outc:
        outc.write(f"# Tailocin Finder Batch Report — {n} samples — {timestamp}\n")
        outc.write("# ANSI-colored version. View with `less -R` or `cat` in a color-capable terminal.\n\n")
        
        outc.write("## SUMMARY\n")
        outc.write("No.\tsample\tcandidates\tbinary\tcandidate_spans\tcluster_sizes_kb\treport_dir\n")
        
        for i, it in enumerate(items, 1):
            spans = it.get("candidate_spans", "")
            sizes = it.get("cluster_sizes_kb", "")
            btxt = f"{it['binary']}"
            btxt_col = f"{ANSI_GREEN}{btxt}{ANSI_RESET}" if it['binary'] == 1 else f"{ANSI_RED}{btxt}{ANSI_RESET}"
            # Optional: If size>10 kb, highlight in green (example)
            sizes_col = sizes if not sizes or float(sizes.split(';')[0]) <= 10 else f"{ANSI_GREEN}{sizes}{ANSI_RESET}"
            outc.write(f"{i}\t{it['sample']}\t{it['n_cand']}\t{btxt_col}\t{spans}\t{sizes_col}\t{it['dir']}\n")
        
        outc.write("\n## DETAILS\n\n")
        
        for i, it in enumerate(items, 1):
            outc.write(f"### [{i}] Sample: {it['sample']}  ({it['dir']})\n\n")
            try:
                with open(it["report_txt"], "r", encoding="utf-8") as rf:
                    content = rf.read()
                colored_content = colorize_pass_fail_block(content)
                outc.write(colored_content.rstrip() + "\n\n")
            except Exception as e:
                outc.write(f"{ANSI_RED}[ERROR]{ANSI_RESET} Cannot read {it['report_txt']}: {e}\n\n")

    return plain_path, colored_path


def collect_reports(results_root="results", out_tsv="batch_report.tsv", out_json="batch_report.json", do_print=False):
    records = []
    for root, dirs, files in os.walk(results_root):
        if "tailocin_report.json" in files:
            jpath = os.path.join(root, "tailocin_report.json")
            try:
                with open(jpath, "r", encoding="utf-8") as f:
                    j = json.load(f)
                params = j.get("params", {})
                
                cand_spans = []
                cand_sizes_kb = []  # New: Used to store the size (kb) of each candidate cluster
                cand_funcs = []     # Keep the original functions
                
                for c in j.get("clusters", []):
                    if c.get("is_candidate"):
                        start = c.get("start")
                        end = c.get("end")
                        if start is not None and end is not None:
                            span_str = f"{start}..{end}"
                            size_kb = (end - start + 1) / 1000.0
                            cand_spans.append(span_str)
                            cand_sizes_kb.append(f"{size_kb:.2f}")
                        cand_funcs.append(",".join(c.get("functions", [])))
                
                rec = {
                    "sample": j.get("sample", os.path.basename(root)),
                    "candidates": j.get("candidates", 0),
                    "binary": j.get("binary_conclusion", 0),
                    "block_scope": j.get("block_scope", ""),
                    "window": params.get("window", ""),
                    "min_cluster_span": params.get("min_cluster_span", ""),
                    "pident": params.get("pident", ""),
                    "length": params.get("length", ""),
                    "bitscore": params.get("bitscore", ""),
                    "evalue": params.get("evalue", ""),
                    "candidate_spans": ";".join(cand_spans),          # New added in 20260303
                    "cluster_sizes_kb": ";".join(cand_sizes_kb),  # New added in 20260303
                    "candidate_functions": ";".join(cand_funcs),
                    "report_dir": root
                }
                records.append(rec)
            except Exception:
                continue
    
    records.sort(key=lambda r: r["sample"])
    
    # Modify headers: Add cluster_Sizes_kb and remove the contig part from the original candidate_stans (doon in 20260303)
    headers = [
        "sample", "candidates", "binary", "block_scope", "window", "min_cluster_span",
        "pident", "length", "bitscore", "evalue", 
        "candidate_spans", "cluster_sizes_kb", "candidate_functions", "report_dir"
    ]
    
    with open(out_tsv, "w", encoding="utf-8") as tf:
        tf.write("\t".join(headers) + "\n")
        for r in records:
            row = [str(r.get(h, "")) for h in headers]
            tf.write("\t".join(row) + "\n")
    
    agg = {
        "total_samples": len(records),
        "positives": sum(1 for r in records if int(r["candidates"]) > 0),
        "negatives": sum(1 for r in records if int(r["candidates"]) == 0),
        "records": records
    }
    with open(out_json, "w", encoding="utf-8") as jf:
        json.dump(agg, jf, indent=2, ensure_ascii=False)
    
    if do_print:
        print(f"Collected {len(records)} sample(s). Positives={agg['positives']}, Negatives={agg['negatives']}")
        for r in records:
            mark = colored("POS", "green") if int(r["candidates"]) > 0 else colored("NEG", "red")
            print(f"[{mark}] {r['sample']}  candidates={r['candidates']}  spans={r['candidate_spans']}  sizes_kb={r['cluster_sizes_kb']}")
    
    # Return records, in order to use aggregate_batch_text
    return records 


# --------------------------- CLI ---------------------------

EXAMPLES_TEXT = """\
Examples:
  # single sample(FASTA), 15kb window, cluster span ≥13.4kb
  python TattleTail.py genome.fna -o results --window 15000 --min-cluster-span 13400

  # Batch process (mix .fna & .ffn)
  python TattleTail.py assemblies/ dir_of_ffn/cds_from_genomic.ffn -o results --window 15000

  # Specify database path and index prefix (relative path, without suffix)
  python TattleTail.py genome.fna -o results --db-fasta database_tailocin.fasta --db-name database_tailocin

  # Print each sample on the terminal (including color) and export the hit details TSV
  python TattleTail.py genome.fna -o results --per-sample-stdout --force-stdout --dump-hits-tsv

  # Re-summary and report (scan results root directory)
  python TattleTail.py collect --results-root results --out-tsv all.tsv --out-json all.json --print
"""

class SmartHelp(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def build_parser():
    p = argparse.ArgumentParser(
        prog="tailocin_finder.py",
        description="Tailocin Finder (cluster-aware) — identify tailocin-encoding BGCs from genomes or FFN.",
        formatter_class=SmartHelp,
        epilog=EXAMPLES_TEXT
    )

    # analysis
    p.add_argument("inputs", nargs="*", help="Input files or directories (.fna/.fa/.fasta or .ffn)")
    p.add_argument("-o", "--outdir", default="results", help="Root output directory")

    # database
    p.add_argument("--db-fasta", default="database_tailocin.fasta", help="Tailocin DB FASTA (for label/function mapping)")
    p.add_argument("--db-name",  default="database_tailocin",       help="BLAST DB prefix/name (blastp -db)")

    # blast threshold
    p.add_argument("--evalue",   type=float, default=1e-10, help="E-value cutoff (also passed to blastp)")
    p.add_argument("--pident",   type=float, default=35.0,  help="Identity cutoff (percent)")
    p.add_argument("--length",   type=int,   default=50,    help="Alignment length cutoff (aa)")
    p.add_argument("--bitscore", type=float, default=50.0,  help="Bitscore cutoff")
    p.add_argument("--threads",  type=int,   default=None,  help="blastp threads")

    # clustering
    p.add_argument("--window", type=int, default=10000, help="Cluster window on same contig (bp)")
    p.add_argument("--min-cluster-span", type=int, default=0, help="Minimum cluster total span to keep (bp). 13400 for strict")

    # logic
    p.add_argument("--block-scope", choices=["cluster","global"], default="cluster", help="Blocker scope: within cluster or genome-wide")
    p.add_argument("--binary-output", action="store_true", help="Use binary wording in conclusion (0/1)")

    # mode
    p.add_argument("--prodigal-mode", choices=["meta","single"], default="meta", help="Prodigal mode for unannotated genomes")
    p.add_argument("--genetic_code", type=int, default=11, help="NCBI translation table (11 for bacteria/archaea)")
    p.add_argument("--ffn-partial", choices=["trim","pad","skip","error"], default="trim",
                   help="When FFN CDS length mod 3 != 0: trim trailing, pad N, skip, or error")

    # markers
    p.add_argument("--flanks",     default="trpE,trpG",                  help="Comma-separated flank markers")
    p.add_argument("--regulators", default="prtN,prtR",                  help="Comma-separated regulator markers")
    p.add_argument("--toxins",     default="holin,endolysin",            help="Comma-separated lysis markers")
    p.add_argument("--blockers",   default="capsid,terminase,integrase", help="Comma-separated phage blockers")

    # outputs
    p.add_argument("--summary-file", default="", help="Append per-sample summary to this TSV (sample, candidates, binary)")
    p.add_argument("--per-sample-stdout", action="store_true", help="Echo per-sample report to stdout (colored)")
    p.add_argument("--force-stdout", action="store_true", help="Force stdout even under non-tty (nohup)")
    p.add_argument("--quiet", action="store_true", help="Do not print reports to stdout")
    p.add_argument("--dump-hits-tsv", action="store_true", help="Write hits_with_coords.tsv for each sample")

    # version
    p.add_argument("--version", action="version", version="TattleTail(VERSION1.0)")

    return p

def parse_args():
    if len(sys.argv) > 1 and sys.argv[1] == "collect":
        cparser = argparse.ArgumentParser(
            prog="tailocin_finder.py collect",
            description="Collect per-sample JSON reports under a results root and write a batch summary (TSV/JSON).",
            formatter_class=SmartHelp,
            epilog=EXAMPLES_TEXT
        )
        cparser.add_argument("--results-root", default="results", help="Root directory to search recursively (use the same as --outdir when running)")
        cparser.add_argument("--out-tsv", default="batch_report.tsv", help="Output TSV path")
        cparser.add_argument("--out-json", default="batch_report.json", help="Output JSON path")
        cparser.add_argument("--print", dest="do_print", action="store_true", help="Print brief summary to stdout")
        
       
        args = cparser.parse_args(sys.argv[2:])
        return args, "collect"

    parser = build_parser()
    # print help without parameters
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    return parser.parse_args(), "run"

# --------------------------- main process ---------------------------

def write_hits_tsv(path, ghits):
    with open(path, "w", encoding="utf-8") as f:
        f.write("qseqid\tcontig\tstart\tend\tfunc\tlabel\tpident\tlength\tbitscore\tevalue\n")
        for h in ghits:
            f.write("{q}\t{c}\t{s}\t{e}\t{func}\t{lab}\t{pi:.3f}\t{al}\t{bs:.1f}\t{ev:.2e}\n".format(
                q=h["q"], c=h.get("contig",""), s=h.get("start",""), e=h.get("end",""),
                func=h["func"], lab=h["label"], pi=h["pident"], al=h["length"],
                bs=h["bitscore"], ev=h["evalue"]
            ))

def write_cluster_annotation_tsv(path, clusters_eval):
    """
    Write detailed annotation table for candidate clusters only.
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write("cluster_id\tcontig\tstart\tend\tgene_id\tfunction\tpident\tlength\tbitscore\tevalue\n")

        for idx, c in enumerate(clusters_eval, 1):
            if not c["is_candidate"]:
                continue

            # Sort by genomic location
            members_sorted = sorted(c["members"], key=lambda x: x["start"])

            for m in members_sorted:
                f.write(
                    f"{idx}\t"
                    f"{c['contig']}\t"
                    f"{m['start']}\t"
                    f"{m['end']}\t"
                    f"{m['q']}\t"
                    f"{m['func']}\t"
                    f"{m['pident']:.2f}\t"
                    f"{m['length']}\t"
                    f"{m['bitscore']:.1f}\t"
                    f"{m['evalue']:.2e}\n"
                )

def analyze_file(
    input_file, output_dir, db_name, db_fasta,
    flank_markers, regulatory_markers, toxin_markers, phage_blockers,
    evalue_thr, pident_thr, length_thr, bitscore_thr, threads,
    prodigal_mode, genetic_code, ffn_partial_policy,
    window_bp, min_cluster_span, block_scope, binary_output,
    quiet_flag=False, force_stdout_flag=False, per_sample_stdout=False, summary_file_path="",
    dump_hits_tsv=False
):
    os.makedirs(output_dir, exist_ok=True)
    sample_name = os.path.basename(os.path.splitext(input_file)[0])
    query_faa  = os.path.join(output_dir, "query_proteins.faa")
    blast_out  = os.path.join(output_dir, "blast_results.txt")
    report_txt = os.path.join(output_dir, "tailocin_report.txt")
    report_json= os.path.join(output_dir, "tailocin_report.json")
    run_log    = os.path.join(output_dir, "run.log")

    with open_log_append(run_log) as lf:
        lf.write(f"[INFO] Analyzing: {input_file}\n")
        lf.write(f"[INFO] Output: {output_dir}\n")
        lf.write(f"[INFO] DB name: {db_name} | DB fasta: {db_fasta}\n")

    annotated = is_annotated_fasta(input_file)
    with open_log_append(run_log) as lf:
        lf.write("[INFO] Annotated FASTA detected; translating CDS directly.\n" if annotated
                 else f"[INFO] Unannotated FASTA; running Prodigal ({prodigal_mode}).\n")

    if annotated:
        try:
            translate_ffn_to_faa(input_file, query_faa, genetic_code=genetic_code,
                                 partial_policy=ffn_partial_policy, log_path=run_log)
        except Exception as e:
            if not quiet_flag:
                print(f"[ERROR] FFN translate failed: {e} (see {run_log})")
            if summary_file_path:
                with open(summary_file_path, "a", encoding="utf-8") as sf:
                    sf.write(f"{sample_name}\t-1\t0\n")
            return {"ok": False, "sample": sample_name, "report_txt": report_txt, "n_cand": -1, "binary": 0, "dir": output_dir}
    else:
        ok, msg = run_prodigal(input_file, query_faa, prodigal_mode=prodigal_mode, genetic_code=genetic_code, log_path=run_log)
        if not ok:
            if not quiet_flag:
                print(f"[ERROR] {msg} (see {run_log})")
            if summary_file_path:
                with open(summary_file_path, "a", encoding="utf-8") as sf:
                    sf.write(f"{sample_name}\t-1\t0\n")
            return {"ok": False, "sample": sample_name, "report_txt": report_txt, "n_cand": -1, "binary": 0, "dir": output_dir}

    ok, msg = run_blastp(query_faa, db_name, blast_out, evalue=evalue_thr, threads=threads, log_path=run_log)
    if not ok:
        if not quiet_flag:
            print(f"[ERROR] {msg} (see {run_log})")
        if summary_file_path:
            with open(summary_file_path, "a", encoding="utf-8") as sf:
                sf.write(f"{sample_name}\t-1\t0\n")
        return {"ok": False, "sample": sample_name, "report_txt": report_txt, "n_cand": -1, "binary": 0, "dir": output_dir}

    id_to_function, id_to_label = parse_functional_mapping(db_fasta)
    hits, found_functions_global = parse_all_hits(
        blast_out, id_to_function, id_to_label,
        evalue_thr=evalue_thr, pident_thr=pident_thr, length_thr=length_thr, bitscore_thr=bitscore_thr
    )

    clusters_eval = []

    # global exclude(Optional)--This branch remains, but simplifies processing
    if block_scope == "global":
        has_global_blocker = any(b.lower() in found_functions_global for b in phage_blockers)
        if has_global_blocker:
            # If there are blockers globally, report directly that there are no candidates
            params = {
                "evalue": evalue_thr, "pident": pident_thr, "length": length_thr, "bitscore": bitscore_thr,
                "window": window_bp, "min_cluster_span": (0 if min_cluster_span <= 0 else min_cluster_span),
                "block_scope": block_scope
            }
            write_report_and_json(
                report_txt, report_json,
                sample_name, params,
                found_functions_global, clusters_eval,  # clusters_eval still null value
                flank_markers, regulatory_markers, toxin_markers, phage_blockers,
                block_scope, binary_output
            )
            # You can choose not to generate tsv because there are no candidates available
            n_cand = 0
            binary = 0
            if summary_file_path:
                with open(summary_file_path, "a", encoding="utf-8") as sf:
                    sf.write(f"{sample_name}\t{n_cand}\t{binary}\n")
            should_print = per_sample_stdout and (sys.stdout.isatty() or force_stdout_flag) and (not quiet_flag)
            if should_print:
                print_report_to_console(sample_name, params, found_functions_global, clusters_eval,
                                        flank_markers, regulatory_markers, toxin_markers, phage_blockers,
                                        block_scope, binary_output)
            return {"ok": True, "sample": sample_name, "report_txt": report_txt, "n_cand": n_cand, "binary": binary, "dir": output_dir}

    # Normal clustering process (regardless of whether block_stope is a cluster or global but without a blocker)
    ghits = []
    if window_bp and window_bp > 0:
        qcoords = parse_prodigal_coords(query_faa)
        if not qcoords:
            qcoords = parse_ffn_header_coords(query_faa, log_path=run_log)
            if qcoords:
                with open_log_append(run_log) as lf:
                    lf.write("[INFO] Coordinates parsed from FFN-style headers.\n")
        if not qcoords:
            with open_log_append(run_log) as lf:
                lf.write("[WARN] No coordinates parsed (Prodigal/FFN); clustering skipped.\n")
        else:
            ghits = join_hits_with_coords(hits, qcoords)
            if dump_hits_tsv and ghits:
                write_hits_tsv(os.path.join(output_dir, "hits_with_coords.tsv"), ghits)
                with open_log_append(run_log) as lf:
                    lf.write(f"[INFO] Dumped hits_with_coords.tsv with {len(ghits)} rows.\n")
            if not ghits:
                with open_log_append(run_log) as lf:
                    lf.write("[WARN] No hits with coordinates after filtering; clustering skipped.\n")
            else:
                clusters = cluster_by_window(ghits, window_bp=window_bp)
                if min_cluster_span and min_cluster_span > 0:
                    before = len(clusters)
                    clusters = [c for c in clusters if (c["end"] - c["start"] + 1) >= min_cluster_span]
                    after = len(clusters)
                    with open_log_append(run_log) as lf:
                        lf.write(f"[INFO] Filtered clusters by min span {min_cluster_span} bp: {before} -> {after}\n")
                clusters_eval = evaluate_clusters_strict(
                    clusters,
                    flank=tuple(flank_markers),
                    blockers=tuple(phage_blockers),
                    regulators=tuple(regulatory_markers),
                    toxins=tuple(toxin_markers)
                )
                # trim flanks
                clusters_eval = trim_cluster_by_flanks(clusters_eval, flank_markers)

    params = {
        "evalue": evalue_thr, "pident": pident_thr, "length": length_thr, "bitscore": bitscore_thr,
        "window": window_bp, "min_cluster_span": (0 if min_cluster_span <= 0 else min_cluster_span),
        "block_scope": block_scope
    }

    # write report
    write_report_and_json(
        report_txt, report_json,
        sample_name, params,
        found_functions_global, clusters_eval,
        flank_markers, regulatory_markers, toxin_markers, phage_blockers,
        block_scope, binary_output
    )

    # New: generate cluster_annotation.tsv(only have)
    n_cand = sum(1 for c in clusters_eval if c["is_candidate"])
    if n_cand > 0:
        annotation_path = os.path.join(output_dir, "cluster_annotation.tsv")
        write_cluster_annotation_tsv(annotation_path, clusters_eval)
        with open_log_append(run_log) as lf:
            lf.write(f"[INFO] Wrote cluster_annotation.tsv with {n_cand} candidate clusters\n")
    else:
        with open_log_append(run_log) as lf:
            lf.write("[INFO] No candidate clusters found, skipping cluster_annotation.tsv\n")

    binary = 1 if n_cand > 0 else 0
    if summary_file_path:
        with open(summary_file_path, "a", encoding="utf-8") as sf:
            sf.write(f"{sample_name}\t{n_cand}\t{binary}\n")

    should_print = per_sample_stdout and (sys.stdout.isatty() or force_stdout_flag) and (not quiet_flag)
    if should_print:
        print_report_to_console(sample_name, params, found_functions_global, clusters_eval,
                                flank_markers, regulatory_markers, toxin_markers, phage_blockers,
                                block_scope, binary_output)
    # calculate spans & sizes
    cand_spans = []
    cand_sizes_kb = []
    for c in clusters_eval:
        if c["is_candidate"]:
            start = c.get("start")
            end = c.get("end")
            if start is not None and end is not None:
                span = f"{start}..{end}"
                size_kb = (end - start + 1) / 1000.0
                cand_spans.append(span)
                cand_sizes_kb.append(f"{size_kb:.2f}")

    return {
        "ok": True,
        "sample": sample_name,
        "report_txt": report_txt,
        "n_cand": n_cand,
        "binary": binary,
        "dir": output_dir,
        "candidate_spans": ";".join(cand_spans),          # 20260303
        "cluster_sizes_kb": ";".join(cand_sizes_kb)       # 20260303
    }





def list_input_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            for fname in os.listdir(p):
                if fname.endswith((".fasta", ".fna", ".ffn", ".fa")):
                    files.append(os.path.join(p, fname))
        elif p.endswith((".fasta", ".fna", ".ffn", ".fa")):
            files.append(p)
        else:
            print(f"[WARNING] Skipping unsupported input: {p}")
    return files

def main():
    args, mode = parse_args()

    if mode == "collect":
        collect_reports(results_root=args.results_root, out_tsv=args.out_tsv, out_json=args.out_json, do_print=args.do_print)
        return

    # Dependency check & BLAST library self-test
    check_dependencies(require_prodigal=True, require_blastp=True)
    ensure_blast_db(db_fasta=args.db_fasta, db_name=args.db_name)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    results_root = args.outdir
    os.makedirs(results_root, exist_ok=True)

    flank_markers      = split_csv(args.flanks)
    regulatory_markers = split_csv(args.regulators)
    toxin_markers      = split_csv(args.toxins)
    phage_blockers     = split_csv(args.blockers)

    processed = []

    input_files = list_input_files(args.inputs)
    if not input_files:
        print("[BATCH] No valid inputs found.")
        sys.exit(1)

    for input_file in input_files:
        sample_name = os.path.splitext(os.path.basename(input_file))[0]
        output_dir = os.path.join(results_root, f"{sample_name}_{timestamp}")
        rec = analyze_file(
            input_file, output_dir, args.db_name, args.db_fasta,
            flank_markers, regulatory_markers, toxin_markers, phage_blockers,
            args.evalue, args.pident, args.length, args.bitscore, args.threads,
            args.prodigal_mode, args.genetic_code, args.ffn_partial,
            args.window, args.min_cluster_span, args.block_scope, args.binary_output,
            quiet_flag=args.quiet, force_stdout_flag=args.force_stdout,
            per_sample_stdout=args.per_sample_stdout, summary_file_path=args.summary_file,
            dump_hits_tsv=args.dump_hits_tsv
        )
        if rec: processed.append(rec)

    if processed:
        plain_path, colored_path = aggregate_batch_text(results_root, processed, timestamp)
        if not args.quiet:
            pos = sum(1 for r in processed if r.get("n_cand", 0) > 0)
            neg = sum(1 for r in processed if r.get("n_cand", 0) == 0)
            print(f"[BATCH] Samples processed: {len(processed)}  Positives: {pos}  Negatives: {neg}")
            print(f"[BATCH] Aggregated report (plain)  → {plain_path}")
            print(f"[BATCH] Aggregated report (colored)→ {colored_path}")
            print("       View colored with: less -R <path>  (or cat)")
    else:
        if not args.quiet:
            print("[BATCH] No valid inputs processed.")

if __name__ == "__main__":
    main()

