# TattleTail: Pyocin Prediction Tool (v1.0)

TattleTail is a command-line pipeline for identifying **canonical candidate pyocin-encoding regions** in *Pseudomonas aeruginosa* genome assemblies. It combines gene prediction, BLASTP searches, genomic coordinates, and rule-based filtering.

Version 1.0 evaluates:

- the flanking genes `trpE` and `trpG`;
- the regulatory genes `prtN` and `prtR`;
- the lysis genes `holin` and `endolysin`; 
- qualifying `capsid`-, `terminase`-, or `integrase`-like blocker hits.


## Scope and interpretation

TattleTail v1.0 is designed to identify canonical candidates occurring between `trpE` and `trpG`. Therefore:

- the required markers must be detected on the same contig;
- non-canonical loci outside this genomic context may not be detected;
- a locus split across contigs may produce a false-negative result;
- a prediction does not demonstrate that a locus is functional; 
- “no candidate found” means that no region passed the selected rules and parameters, not necessarily that the genome lacks all pyocin-related genes.

## Workflow

For each genome, TattleTail:

1. predicts proteins with Prodigal;
2. searches them against the TattleTail protein database using BLASTP;
3. filters hits by identity, alignment length, E-value and bit score;
4. joins qualifying hits to genomic coordinates;
5. groups consecutive hits on the same contig when their separation does not exceed `--window`;
6. evaluates each preliminary hit cluster for a valid `trpE` and `trpG` flanking pair and defines the inter-flank core when such a pair is available;
7. evaluates the preliminary hit cluster for `prtN`, `prtR`, `holin`, and `endolysin`;
8. applies the selected phage-blocker rule to the preliminary hit cluster;
9. applies `--min-cluster-span` to the trimmed inter-flank core when that core can be defined.

Under the default cluster-level rule, one qualifying `capsid`-, `terminase`-, or `integrase`-like BLASTP hit within a preliminary hit cluster is sufficient for that cluster to fail the phage-blocker criterion, regardless of whether an inter-flank core can be defined.

## Requirements

- Python 3.10.20
- Biopython 1.87
- termcolor 3.3.0
- Prodigal 2.6.3
- NCBI BLAST+ 2.17.0 
- Bakta 1.12.0 and a compatible Bakta database (light v6.0 - 2025-02-24)

Other compatible versions may also work but have not been tested systematically.

## Installation

### Clone the repository

```bash
git clone https://github.com/xthua/tattletail.git
cd tattletail
```

### Create a reproducible Conda environment

The repository provides an `environment.yml` file containing the tested dependency versions. Create the complete software environment with one command:

```bash
conda env create -f environment.yml
conda activate tattletail_env
```

Alternatively, the same environment can be created directly from the command line:

```bash
conda create -n tattletail_env \
    -c conda-forge -c bioconda \
    --strict-channel-priority \
    python=3.10.20 \
    biopython=1.87 \
    termcolor=3.3.0 \
    prodigal=2.6.3 \
    blast=2.17.0 \
    bakta=1.12.0

conda activate tattletail_env
```

Check the installation:

```bash
python --version
prodigal -v
blastp -version
makeblastdb -version
bakta --version
```

On Linux, the following command can be used to confirm that the executables are being obtained from the activated Conda environment:

```bash
which python prodigal blastp makeblastdb bakta
```

Installing Bakta does not automatically install its annotation database.

## Database preparation

### TattleTail protein database

The repository must contain the supplied `database_tailocin.fasta`. Build and verify the BLAST database:

```bash
makeblastdb \
    -in database_tailocin.fasta \
    -dbtype prot \
    -out database_tailocin

blastdbcmd -db database_tailocin -info
```

### Bakta database

```bash
bakta_db download --output bakta_db_light --type light
```

The examples below assume that the database is at `bakta_db_light/db-light`. Otherwise, supply its path with `--bakta-db`.

See the [official Bakta database instructions](https://github.com/oschwengers/bakta#database) for database installation and version compatibility.

## Reproducible PAO1 test

This test uses the complete *P. aeruginosa* PAO1 reference assembly.

| Item | Expected value |
|---|---|
| Assembly | GCF_000006765.1 |
| Input | GCF_000006765.1_ASM676v1_genomic.fna |
| candidate count | 1 |
| candidate_spans | 672459..703476 |
| cluster_sizes_kb | 31.02 kb |

Download the exact input:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz

gunzip GCF_000006765.1_ASM676v1_genomic.fna.gz
```

Run the test:

```bash
python TattleTail.py \
    GCF_000006765.1_ASM676v1_genomic.fna \
    -o PAO1_test \
    --db-name database_tailocin \
    --db-fasta database_tailocin.fasta \
    --bakta-db bakta_db_light/db-light \
    --window 15000 \
    --min-cluster-span 13400 
```

This smoke test reproduces the settings used in the reported analyses. The program default for `--min-cluster-span` is `0`.

The summary should report:

```text
candidates: 1
candidate_spans: 672459..703476
cluster_sizes_kb: 31.02
```

Coordinates are 1-based and inclusive and describe the final trimmed inter-flank core.

```bash
less -R PAO1_test/GCF_000006765.1_ASM676v1_genomic_*/tailocin_report.txt
```

If the result differs, confirm the exact assembly, input file, TattleTail version, reference database, Bakta database, and command. Different PAO1 records or database versions may produce different results.

## Usage

Single or multiple genomes:

```bash
python TattleTail.py genome.fna -o results
python TattleTail.py genome1.fna genome2.fna -o results
```

Batch input:

```bash
python TattleTail.py dataset/ -o results
python TattleTail.py dataset/*.fna -o results
python TattleTail.py dataset/*/*.fna -o results
```

Mixed genomic FASTA and nucleotide CDS FASTA (FFN) input:

```bash
python TattleTail.py assemblies/*.fna dataset_ffn/*/*.ffn -o results
```

Background run:

```bash
nohup python TattleTail.py dataset/*/*.fna \
    -o results \
    --window 15000 \
    --min-cluster-span 13400 \
    > tattletail.nohup.log 2>&1 &
```

## Parameters

| Option | Program default | Reported analyses | Meaning |
|---|---:|---:|---|
| `--evalue` | `1e-10` | `1e-10` | Maximum BLASTP E-value |
| `--pident` | `35` | `35` | Minimum amino-acid identity (%) |
| `--length` | `50` | `50` | Minimum alignment length (aa) |
| `--bitscore` | `50` | `50` | Minimum BLASTP bit score |
| `--window` | `15000` | `15000` | Maximum separation between consecutive qualifying hits |
| `--min-cluster-span` | `0` | `13400` | Minimum trimmed core length; `0` disables the filter |
| `--block-scope` | `cluster` | `cluster` | Scope of the phage-blocker rule |
| `--prodigal-mode` | `meta` | `meta` | Prodigal mode |
| `--genetic_code` | `11` | `11` | Prodigal translation table |
| `--bakta-db` | `bakta_db_light/db-light` | database-dependent | Bakta database path |

Additional input, marker and output-control options:

| Option | Default | Meaning |
|---|---:|---|
| `--threads` | automatic | Number of BLASTP threads |
| `--binary-output` | disabled | Use binary wording in the report conclusion |
| `--ffn-partial` | `trim` | Handling of partial codons in FFN input: `trim`, `pad`, `skip`, or `error` |
| `--flanks` | `trpE,trpG` | Comma-separated flank-marker labels |
| `--regulators` | `prtN,prtR` | Comma-separated regulator-marker labels |
| `--toxins` | `holin,endolysin` | Comma-separated lysis-marker labels |
| `--blockers` | `capsid,terminase,integrase` | Comma-separated phage-blocker labels |
| `--summary-file` | none | Append the per-sample summary to a specified TSV file |
| `--per-sample-stdout` | disabled | Print each sample report to standard output |
| `--force-stdout` | disabled | Print reports even when standard output is not a terminal |
| `--quiet` | disabled | Suppress report printing to standard output |
| `--dump-hits-tsv` | disabled | Write `hits_with_coords.tsv` for qualifying hits |
| `--version` | — | Display the TattleTail version and exit |

```bash
python TattleTail.py --help
```

### Window

The 15,000-bp default was selected empirically using 98 genomes with manually identified canonical pyocin regions. A 10,000-bp window fragmented marker-hit chains, whereas 15,000 bp retained the complete regions. Larger windows can connect more distant markers but can also merge unrelated phage-derived regions. Regions recovered only with a larger window require cautious manual interpretation.

### Minimum cluster span

For the reported analyses, 13,400 bp was set below the smallest positive core observed among the 98 positive genomes analysed with a 15,000-bp window. It suppresses clearly short, fragmented, or redundant clusters. This is an empirical reporting filter, not a universal biological minimum. The program default is `0`.

### Blocker scope

- `--block-scope cluster`: checks each preliminary hit cluster. One qualifying blocker causes that cluster to fail the phage-blocker criterion.
- `--block-scope global`: checks qualifying blocker hits across the entire genome. A qualifying blocker anywhere in the genome prevents candidate reporting under this mode.

## Input guidance

- Genomic nucleotide FASTA files (`.fna`, `.fa`, or `.fasta`) are recommended.
- Nucleotide CDS FASTA files (`.ffn`) are supported when their headers contain usable sequence and coordinate information.
- Complete genomes are preferred for reproducible coordinate-level results.
- Draft assemblies can be analysed, but candidates cannot be reconstructed across contigs.
- Prefer simple, unique FASTA identifiers.
- FFN input can be useful for pre-annotated data, but it is not equivalent to a contiguous genome assembly for evaluating locus continuity.

## Outputs

Each sample receives a timestamped directory.

| File or directory | Description |
|---|---|
| `tailocin_report.txt` | Human-readable report |
| `tailocin_report.json` | Machine-readable report |
| `blast_results.txt` | Raw BLASTP output |
| `query_proteins.faa` | Prodigal proteins |
| `run.log` | Commands, parameters, and run information |
| `bakta_full_genome/` | Bakta annotation |

Accepted candidates can also produce:

| File | Description |
|---|---|
| `cluster_N_gene_list.tsv` | Genes in candidate N |
| `cluster_N_on_full_contig.gff3` | Candidate features in GFF3 |
| `cluster_N_on_full_contig.gbk` | Candidate features in GenBank format |
| `cluster_annotation.tsv` | Consolidated annotation |

Batch summaries:

- `allresultsof<N>samples_<timestamp>.txt`
- `allresultsof<N>samples_<timestamp>_colored.txt`

Key fields are `sample`, `candidates`, `binary`, `candidate_spans`, and `cluster_sizes_kb`. `cluster_sizes_kb` is calculated from inclusive coordinates and rounded to two decimal places.

## Collect existing results

```bash
python TattleTail.py collect \
    --results-root results \
    --out-tsv all_results.tsv \
    --out-json all_results.json \
    --print
```

This rebuilds a batch summary without rerunning genomes.

If `--out-tsv` and `--out-json` are omitted, the default output names are `batch_report.tsv` and `batch_report.json`, respectively.

## Export qualifying hits

```bash
python TattleTail.py genome.fna -o results --dump-hits-tsv
```

This table helps explain cluster breaks and blocker decisions.

## Troubleshooting

| Problem | Checks and possible solutions |
|---|---|
| `prodigal: command not found` | Activate the Conda environment and confirm that `which prodigal` points inside that environment. |
| `blastp` or `makeblastdb: command not found` | Confirm that NCBI BLAST+ is installed in the activated environment. |
| `bakta: command not found` | Activate the Conda environment and run `bakta --version`. |
| Bakta database not found | Download the database and provide the extracted `db-light` directory with `--bakta-db`. |
| Bakta software/database incompatibility | Compare `bakta --version` with the Bakta database information in `version.json`. |
| BLAST database files are missing or invalid | Rebuild the database using `makeblastdb -in database_tailocin.fasta -dbtype prot -out database_tailocin`. |
| PAO1 result differs | Confirm the exact assembly, versions, databases, and full command. |
| Genome-wide markers but no candidate | Use `--dump-hits-tsv`; a gap larger than `--window` may split the markers. |
| `Unable to define an inter-flank core region` | The cluster lacked a usable `trpE`–`trpG` pair defining a positive interval. |
| Candidate rejected by blockers | Inspect blocker hits within the selected scope. |
| Draft genome produces no candidate | Check whether the candidate region or required flanks are split across contigs; TattleTail does not combine evidence across contig boundaries. |
| Bakta rejects a FASTA header | Use short identifiers without INSDC-incompatible punctuation and retain a mapping file. If necessary, remove `--keep-contig-headers` from the Bakta command in the script. |

## Citation

Rayhaan G. Pais, Weilian Chen, Sebastian Leptihn, Xiaoting Hua, Belinda Loh  
*TattleTail: A Pyocin Prediction Tool.*  
bioRxiv (2026).  
https://doi.org/10.64898/2026.03.25.712926
