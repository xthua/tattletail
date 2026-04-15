#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TattleTail: Pyocin / Tailocin Prediction Tool (v1.0)
"""

import os
import re
import sys
import argparse
import subprocess
import json
import shutil
from shutil import which
from datetime import datetime
from collections import defaultdict

from termcolor import colored
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

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

def check_bakta(bakta_db):
    bakta_path = shutil.which("bakta")
    if bakta_path is None:
        raise RuntimeError(
            "ERROR: Bakta not found in PATH.\n"
            "Please install Bakta (conda install -c conda-forge -c bioconda bakta)."
        )
    if bakta_db is None:
        raise RuntimeError("ERROR: --bakta-db is required.")
    if not os.path.exists(bakta_db) or not os.path.isdir(bakta_db):
        raise RuntimeError(f"ERROR: Bakta database directory not found: {bakta_db}, please download Bakta database by command line.")
    try:
        subprocess.run(
            ["bakta", "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
    except Exception:
        raise RuntimeError(
            "ERROR: Bakta detected but failed to run.\n"
            "Please check your installation."
        )

def parse_gff3_genes(gff_file):
    genes = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = {}
            for x in cols[8].split(";"):
                if "=" in x:
                    k, v = x.split("=", 1)
                    attrs[k] = v
            gene = attrs.get("gene", attrs.get("locus_tag", attrs.get("ID", "")))
            product = attrs.get("product", "Hypothetical protein")
            genes.append({
                "gene": gene,
                "product": product,
                "start": start,
                "end": end,
                "strand": strand
            })
    return genes

def write_cluster_gene_list(outfile, genes):
    """(2) requirement: all genes table (exact format)"""
    with open(outfile, "w") as f:
        f.write("#\tGene\tGene product\tGenomic position\n")
        for i, g in enumerate(genes, 1):
            gene = g.get("gene", "")
            product = g.get("product", "Hypothetical protein")
            start = int(g["start"])
            end = int(g["end"])
            strand = g["strand"]
            start_fmt = f"{start:,}"
            end_fmt = f"{end:,}"
            if strand == "-":
                pos = f"Complement {start_fmt} - {end_fmt}"
            else:
                pos = f"{start_fmt} - {end_fmt}"
            f.write(f"{i}\t{gene}\t{product}\t{pos}\n")

def write_cluster_gff3(outfile, contig_id, genes):
    """Generate standard GFF3 file (best for IGV/JBrowse/Artemis)"""
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("##gff-version 3\n")
        f.write(f"# Tailocin Finder v1.0 - Core pyocin cluster annotations\n")
        f.write(f"# Contig: {contig_id}\n")
        f.write(f"# Only core cluster genes annotated (absolute genomic coordinates)\n\n")
        
        for g in genes:
            start = g["start"]
            end   = g["end"]
            strand = g["strand"]
            gene   = g.get("gene", "unknown")
            product = g.get("product", "Hypothetical protein")
            
            # sanitize attributes to prevent parsing errors
            product_safe = product.replace(';', ',').replace('=', '_').replace('"', "'")
            
            attributes = (
                f"ID={gene};"
                f"Name={gene};"
                f"gene={gene};"
                f"locus_tag={gene};"
                f"product={product_safe};"
                f"Note=Core pyocin cluster gene (Tailocin Finder)"
            )
            
            f.write(f"{contig_id}\tTailocinFinder\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n")

# --------------------------- BLAST DB self-check ---------------------------

def ensure_blast_db(db_fasta="database_tailocin.fasta", db_name="database_tailocin"):
    if not os.path.exists(db_fasta):
        if os.path.basename(db_fasta) == "database_tailocin.fasta":
            print("database_tailocin.fasta does not exist", file=sys.stderr)
        else:
            print(f"[FATAL] DB FASTA does not exist: {db_fasta}", file=sys.stderr)
        sys.exit(2)
    idx_files = [db_name + ".pin", db_name + ".phr", db_name + ".psq"]
    have_index = all(os.path.exists(p) for p in idx_files)
    if not have_index:
        if which("makeblastdb") is None:
            print("[FATAL] makeblastdb not found in PATH", file=sys.stderr)
            print("        Please install BLAST+ (example: conda install -c bioconda blast)", file=sys.stderr)
        cmd = f'makeblastdb -in "{db_fasta}" -dbtype prot -out "{db_name}" -parse_seqids'
        print("[FATAL] BLAST database index not found:", db_name, file=sys.stderr)
        print("        Please create the index first (copy and run the command below):", file=sys.stderr)
        print("        " + cmd, file=sys.stderr)
        sys.exit(2)
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
        print("[FATAL] BLAST database index exists but is unreadable:", db_name, file=sys.stderr)
        print("        Rebuild the index with:", file=sys.stderr)
        print("        " + cmd, file=sys.stderr)
        sys.exit(2)
    try:
        fasta_mtime = os.path.getmtime(db_fasta)
        idx_mtime = min(os.path.getmtime(p) for p in idx_files)
        if fasta_mtime > idx_mtime:
            cmd = f'makeblastdb -in "{db_fasta}" -dbtype prot -out "{db_name}" -parse_seqids'
            print("[WARN] DB FASTA is newer than BLAST index. Rebuild the index to keep consistency:", file=sys.stderr)
            print("       " + cmd, file=sys.stderr)
    except Exception:
        pass
    bad = 0
    try:
        with open(db_fasta, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.strip().split()
                    if len(parts) < 2:
                        bad += 1
        if bad > 0:
            print(f"[FATAL] {db_fasta} contains {bad} headers without a second token (functional label).", file=sys.stderr)
            print("        Headers must follow format: '>ID LABEL [optional]'", file=sys.stderr)
            print("        Example: '>BAR70105.1 integrase_1 [Pseudomonas aeruginosa]'", file=sys.stderr)
            sys.exit(2)
    except Exception:
        pass

# --------------------------- other functions ---------------------------

def check_dependencies(require_prodigal=True, require_blastp=True):
    ok = True
    if require_prodigal and which("prodigal") is None:
        print("prodigal not found in PATH", file=sys.stderr); ok = False
    if require_blastp and which("blastp") is None:
        print("blastp not found in PATH", file=sys.stderr); ok = False
    if not ok:
        sys.exit(1)

def is_annotated_fasta(file_path):
    annotated_keywords = ["gene=", "cds", "product=", "locus_tag", "[location"]
    with open(file_path, 'r', errors="ignore") as f:
        data = f.read(1024)
    return any(k in data.lower() for k in annotated_keywords)

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

def parse_functional_mapping(database_fasta):
    id_to_function = {}
    id_to_label = {}
    with open(database_fasta, errors="ignore") as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:].strip()         
                space_parts = header.split(maxsplit=1)    
                prot_id = space_parts[0]                 
                
                if len(space_parts) > 1 and '~~~' in space_parts[1]:
                    tilde_parts = space_parts[1].split('~~~')
                    label = tilde_parts[0].strip() if tilde_parts else prot_id
                else:
                    label = prot_id
                
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
    def log(msg):
        if log_path:
            with open(log_path, "a") as lf:
                lf.write("[FFN_COORDS] " + msg + "\n")
    qcoords = {}
    for rec in SeqIO.parse(faa_file, "fasta"):
        hdr = rec.description
        token = hdr.split()[0]
        m = re.search(r'^[^|]+\|(.+?)_cds', token)
        if m:
            contig = m.group(1)
        else:
            m2 = re.search(r'^[^|]+\|([^\s]+)', token)
            contig = m2.group(1) if m2 else re.split(r'_cds', token)[0]
        m = re.search(r'\[location=([^\]]+)\]', hdr)
        if not m:
            continue
        loc = m.group(1)
        loc_clean = loc
        for kw in ['complement', 'join', 'order']:
            loc_clean = re.sub(rf'{kw}\(', '(', loc_clean)
        loc_clean = loc_clean.replace('(', '').replace(')', '')
        starts, ends = [], []
        for part in loc_clean.split(','):
            part = part.strip()
            if not part:
                continue
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
    return qcoords

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
    flank_set = set(x.lower() for x in flank_markers)
    for c in cluster_eval:
        if not c["is_candidate"]:
            continue
        flank_hits = [m for m in c["members"] if m["func"].lower() in flank_set]
        if len(flank_hits) < 2:
            continue
        trpE_hit = next((m for m in flank_hits if m["func"].lower() == "trpe"), None)
        trpG_hit = next((m for m in flank_hits if m["func"].lower() == "trpg"), None)
        if not trpE_hit or not trpG_hit:
            continue
        flanks_sorted = sorted([trpE_hit, trpG_hit], key=lambda x: x["start"])
        left_flank = flanks_sorted[0]
        right_flank = flanks_sorted[1]
        new_start = left_flank["end"] + 1
        new_end = right_flank["start"] - 1
        if new_start > new_end:
            continue
        c["start"] = new_start
        c["end"] = new_end
        new_members = []
        for m in c["members"]:
            if m["func"].lower() in flank_set:
                continue
            if m["start"] >= new_start and m["end"] <= new_end:
                new_members.append(m)
        c["members"] = new_members
        c["funcs"] = set(m["func"].lower() for m in new_members)
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

def aggregate_batch_text(results_root, items, timestamp):
    n = len(items)
    base = f"allresultsof{n}samples_{timestamp}"
    plain_path = os.path.join(results_root, base + ".txt")
    colored_path = os.path.join(results_root, base + "_colored.txt")
    os.makedirs(results_root, exist_ok=True)
    with open(plain_path, "w", encoding="utf-8") as out:
        out.write(f"# Tailocin Finder Batch Report — {n} samples — {timestamp}\n")
        out.write("# Each section below is the full content of per-sample tailocin_report.txt.\n\n")
        out.write("## SUMMARY\n")
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
                cand_sizes_kb = []
                cand_funcs = []
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
                    "candidate_spans": ";".join(cand_spans),
                    "cluster_sizes_kb": ";".join(cand_sizes_kb),
                    "candidate_functions": ";".join(cand_funcs),
                    "report_dir": root
                }
                records.append(rec)
            except Exception:
                continue
    records.sort(key=lambda r: r["sample"])
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
    return records

# --------------------------- CLI ---------------------------

EXAMPLES_TEXT = """\
Examples:
  # Single sample (unannotated FASTA): 15kb window, cluster span >=13.4kb, bakta-db: bakta_db_light/db-light(default, please download)
  python TattleTail.py genome.fna -o results --window 15000 --min-cluster-span 13400 --bakta-db db_dir
  # Batch processing (mixed .fna and .ffn)
  python TattleTail.py assemblies/ dataset_ffn/GC*/cds_from_genomic.ffn -o results --window 15000
  # Specify custom database
  python TattleTail.py genome.fna -o results --db-fasta database_tailocin.fasta --db-name database_tailocin
  # Print per-sample reports + export hit details
  python TattleTail.py genome.fna -o results --per-sample-stdout --force-stdout --dump-hits-tsv
  # Re-summarize previous results (no re-run)
  python TattleTail.py collect --results-root results --out-tsv all_results.tsv --out-json all_results.json --print
  # Show all options
  python TattleTail.py -h
"""

class SmartHelp(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

def build_parser():
    p = argparse.ArgumentParser(
        prog="TattleTail.py",
        description="TattleTail — Pyocin / Tailocin Prediction Tool",
        formatter_class=SmartHelp,
        epilog=EXAMPLES_TEXT
    )
    p.add_argument("inputs", nargs="*", help="Input files or directories (.fna/.fa/.fasta or .ffn)")
    p.add_argument("-o", "--outdir", default="results", help="Root output directory")
    p.add_argument("--db-fasta", default="database_tailocin.fasta", help="Tailocin DB FASTA")
    p.add_argument("--db-name",  default="database_tailocin",       help="BLAST DB prefix")
    p.add_argument("--evalue",   type=float, default=1e-10, help="E-value cutoff")
    p.add_argument("--pident",   type=float, default=35.0,  help="Identity cutoff")
    p.add_argument("--length",   type=int,   default=50,    help="Alignment length cutoff")
    p.add_argument("--bitscore", type=float, default=50.0,  help="Bitscore cutoff")
    p.add_argument("--threads",  type=int,   default=None,  help="blastp threads")
    p.add_argument("--window", type=int, default=10000, help="Cluster window (bp)")
    p.add_argument("--min-cluster-span", type=int, default=0, help="Minimum cluster total span (bp)")
    p.add_argument("--block-scope", choices=["cluster","global"], default="cluster", help="Blocker scope")
    p.add_argument("--binary-output", action="store_true", help="Use binary wording in conclusion")
    p.add_argument("--prodigal-mode", choices=["meta","single"], default="meta", help="Prodigal mode")
    p.add_argument("--genetic_code", type=int, default=11, help="NCBI translation table")
    p.add_argument("--ffn-partial", choices=["trim","pad","skip","error"], default="trim", help="FFN partial policy")
    p.add_argument("--flanks",     default="trpE,trpG",                  help="Comma-separated flank markers")
    p.add_argument("--regulators", default="prtN,prtR",                  help="Comma-separated regulator markers")
    p.add_argument("--toxins",     default="holin,endolysin",            help="Comma-separated lysis markers")
    p.add_argument("--blockers",   default="capsid,terminase,integrase", help="Comma-separated phage blockers")
    p.add_argument("--bakta-db", default=os.path.expanduser("bakta_db_light/db-light"), help="Path to Bakta database directory (required)")
    p.add_argument("--summary-file", default="", help="Append per-sample summary to this TSV")
    p.add_argument("--per-sample-stdout", action="store_true", help="Echo per-sample report to stdout (colored)")
    p.add_argument("--force-stdout", action="store_true", help="Force stdout even under non-tty")
    p.add_argument("--quiet", action="store_true", help="Do not print reports to stdout")
    p.add_argument("--dump-hits-tsv", action="store_true", help="Write hits_with_coords.tsv")
    p.add_argument("--version", action="version", version="TattleTail v1.0")
    return p

def parse_args():
    if len(sys.argv) > 1 and sys.argv[1] == "collect":
        cparser = argparse.ArgumentParser(prog="TattleTail.py collect", formatter_class=SmartHelp, epilog=EXAMPLES_TEXT)
        cparser.add_argument("--results-root", default="results", help="Root directory")
        cparser.add_argument("--out-tsv", default="batch_report.tsv", help="Output TSV")
        cparser.add_argument("--out-json", default="batch_report.json", help="Output JSON")
        cparser.add_argument("--print", dest="do_print", action="store_true", help="Print brief summary")
        args = cparser.parse_args(sys.argv[2:])
        return args, "collect"
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args(), "run"

# --------------------------- helper functions ---------------------------

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
    """(3) requirement table: locus tag + alignment statistics per hit"""
    with open(path, "w", encoding="utf-8") as f:
        f.write("cluster_id\tcontig\tstart\tend\tgene_id\tfunction\tpident\tlength\tbitscore\tevalue\n")
        for idx, c in enumerate(clusters_eval, 1):
            if not c["is_candidate"]:
                continue
            members_sorted = sorted(c["members"], key=lambda x: x["start"])
            for m in members_sorted:
                f.write(
                    f"{idx}\t{c['contig']}\t{m['start']}\t{m['end']}\t{m['q']}\t{m['func']}\t"
                    f"{m['pident']:.2f}\t{m['length']}\t{m['bitscore']:.1f}\t{m['evalue']:.2e}\n"
                )

# --------------------------- main analysis function ---------------------------

def analyze_file(
    input_file, output_dir, db_name, db_fasta, bakta_db,
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
    if block_scope == "global":
        has_global_blocker = any(b.lower() in found_functions_global for b in phage_blockers)
        if has_global_blocker:
            params = {
                "evalue": evalue_thr, "pident": pident_thr, "length": length_thr, "bitscore": bitscore_thr,
                "window": window_bp, "min_cluster_span": (0 if min_cluster_span <= 0 else min_cluster_span),
                "block_scope": block_scope
            }
            write_report_and_json(
                report_txt, report_json, sample_name, params,
                found_functions_global, clusters_eval,
                flank_markers, regulatory_markers, toxin_markers, phage_blockers,
                block_scope, binary_output
            )
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

    ghits = []
    if window_bp and window_bp > 0:
        qcoords = parse_prodigal_coords(query_faa)
        if not qcoords:
            qcoords = parse_ffn_header_coords(query_faa, log_path=run_log)
        if qcoords:
            ghits = join_hits_with_coords(hits, qcoords)
            if dump_hits_tsv and ghits:
                write_hits_tsv(os.path.join(output_dir, "hits_with_coords.tsv"), ghits)
            if ghits:
                clusters = cluster_by_window(ghits, window_bp=window_bp)
                if min_cluster_span and min_cluster_span > 0:
                    clusters = [c for c in clusters if (c["end"] - c["start"] + 1) >= min_cluster_span]
                clusters_eval = evaluate_clusters_strict(
                    clusters, flank=tuple(flank_markers), blockers=tuple(phage_blockers),
                    regulators=tuple(regulatory_markers), toxins=tuple(toxin_markers)
                )
                clusters_eval = trim_cluster_by_flanks(clusters_eval, flank_markers)

    # === Full-genome Bakta annotation (run only once per sample) ===
    # This replaces the old per-cluster Bakta calls to avoid boundary drift
    bakta_dir = os.path.join(output_dir, "bakta_full_genome")
    full_gff = os.path.join(bakta_dir, f"{sample_name}.gff3")

    if not os.path.exists(full_gff):
        with open(run_log, "a") as lf:
            lf.write(f"[INFO] Running Bakta on full genome: {input_file}\n")
        cmd = [
            "bakta",
            "--db", bakta_db,
            "--proteins", db_fasta,
            "--output", bakta_dir,
            "--prefix", sample_name,
            "--force",
            "--compliant",
            "--skip-plot",
            "--min-contig-length", "1",
            "--keep-contig-headers",
            input_file
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("=== Bakta ERROR ===")
            print(result.stderr)
            print("=== End of Bakta ERROR ===")
            raise subprocess.CalledProcessError(result.returncode, cmd)

    full_genes = parse_gff3_genes(full_gff)
    # === Load full genome for GBK generation ===
    genome_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))


    # === Filter core genes from full-genome annotation ===
    # 20260410: This replaces the old per-cluster Bakta logic
    for idx, c in enumerate([c for c in clusters_eval if c["is_candidate"]], 1):
        core_start = c["start"]
        core_end   = c["end"]
        
        # Filter genes that fall inside the core cluster
        core_genes = [g for g in full_genes 
                      if core_start <= g["start"] <= core_end 
                      and core_start <= g["end"] <= core_end]
        
        # (2) Write gene_list.tsv
        gene_list_file = os.path.join(output_dir, f"cluster_{idx}_gene_list.tsv")
        write_cluster_gene_list(gene_list_file, core_genes)
        
        # Write full-contig GFF3/GBK (only core genes)
        write_cluster_gff3(
            os.path.join(output_dir, f"cluster_{idx}_on_full_contig.gff3"),
            c["contig"],
            core_genes
        )

        # 20260410: produce cluster_X_on_full_contig.gbk ===
        if c["contig"] in genome_dict:
            full_record = genome_dict[c["contig"]]
            full_record.annotations["molecule_type"] = "DNA"
            full_record.features = []

            for g in core_genes:
                start = g["start"] - 1
                end   = g["end"]
                strand_val = 1 if g["strand"] == "+" else -1
                
                feature = SeqFeature(
                    FeatureLocation(start, end, strand=strand_val),
                    type="CDS",
                    qualifiers={
                        "gene": g.get("gene", "unknown"),
                        "product": g.get("product", "Hypothetical protein"),
                        "locus_tag": g.get("gene", "unknown")
                    }
                )
                full_record.features.append(feature)

            full_gbk = os.path.join(output_dir, f"cluster_{idx}_on_full_contig.gbk")
            SeqIO.write(full_record, full_gbk, "genbank")
    params = {
        "evalue": evalue_thr, "pident": pident_thr, "length": length_thr, "bitscore": bitscore_thr,
        "window": window_bp, "min_cluster_span": (0 if min_cluster_span <= 0 else min_cluster_span),
        "block_scope": block_scope
    }
    write_report_and_json(
        report_txt, report_json, sample_name, params,
        found_functions_global, clusters_eval,
        flank_markers, regulatory_markers, toxin_markers, phage_blockers,
        block_scope, binary_output
    )

    n_cand = sum(1 for c in clusters_eval if c["is_candidate"])

    # (3) BLAST hits statistics table
    if n_cand > 0:
        annotation_path = os.path.join(output_dir, "cluster_annotation.tsv")
        write_cluster_annotation_tsv(annotation_path, clusters_eval)
        with open_log_append(run_log) as lf:
            lf.write(f"[INFO] (3) Wrote cluster_annotation.tsv with {n_cand} candidate clusters\n")

    binary = 1 if n_cand > 0 else 0

    if summary_file_path:
        with open(summary_file_path, "a", encoding="utf-8") as sf:
            sf.write(f"{sample_name}\t{n_cand}\t{binary}\n")

    should_print = per_sample_stdout and (sys.stdout.isatty() or force_stdout_flag) and (not quiet_flag)
    if should_print:
        print_report_to_console(
            sample_name, params, found_functions_global, clusters_eval,
            flank_markers, regulatory_markers, toxin_markers, phage_blockers,
            block_scope, binary_output
        )

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
        "candidate_spans": ";".join(cand_spans),
        "cluster_sizes_kb": ";".join(cand_sizes_kb)
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
    check_bakta(args.bakta_db)
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
            input_file, output_dir, args.db_name, args.db_fasta, args.bakta_db,
            flank_markers, regulatory_markers, toxin_markers, phage_blockers,
            args.evalue, args.pident, args.length, args.bitscore, args.threads,
            args.prodigal_mode, args.genetic_code, args.ffn_partial,
            args.window, args.min_cluster_span, args.block_scope, args.binary_output,
            quiet_flag=args.quiet, force_stdout_flag=args.force_stdout,
            per_sample_stdout=args.per_sample_stdout, summary_file_path=args.summary_file,
            dump_hits_tsv=args.dump_hits_tsv
        )
        if rec:
            processed.append(rec)
    if processed:
        plain_path, colored_path = aggregate_batch_text(results_root, processed, timestamp)
        if not args.quiet:
            pos = sum(1 for r in processed if r.get("n_cand", 0) > 0)
            neg = sum(1 for r in processed if r.get("n_cand", 0) == 0)
            print(f"[BATCH] Samples processed: {len(processed)}  Positives: {pos}  Negatives: {neg}")
            print(f"[BATCH] Aggregated report (plain)  → {plain_path}")
            print(f"[BATCH] Aggregated report (colored)→ {colored_path}")
            print("       View colored with: less -R <path>")
    else:
        if not args.quiet:
            print("[BATCH] No valid inputs processed.")

if __name__ == "__main__":
    main()
