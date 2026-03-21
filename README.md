# TattleTail: Pyocin Prediction Tool (v1.0)

A Python-based bioinformatics tool for identifying potential **tailocin/pyocin** gene clusters in bacterial genomes.

## Features

- Automatic BLAST database self-check and validation
- Supports mixed input formats (FFN or unannotated FASTA)
- Gene prediction using Prodigal or direct translation from FFN
- BLASTP search with functional label mapping
- Window-based clustering of hits
- Strict four-step candidate evaluation
- Detailed per-cluster annotation output
- Batch processing with plain text + colored ANSI summary reports

## Requirements

- **Python** ≥ 3.9
- **Biopython** ≥ 1.85
- **termcolor** ≥ 2.1.0
- **Prodigal** ≥ 2.6.3
- **BLAST+** ≥ 2.16.0 (including `blastp`, `makeblastdb`, `blastdbcmd`)

### Installation

Create and activate a conda environment:

```bash
conda create -n tailocin python=3.9
conda activate tailocin
conda install -c bioconda prodigal blast biopython
pip install termcolor
(If termcolor is already installed via conda, skip the pip step.)

Usage Examples:
Single sample (FASTA input)
15 kb window, minimum cluster span ≥ 13.4 kb
python TattleTail.py genome.fna -o results --window 15000 --min-cluster-span 13400

Batch processing (mixed .fna and .ffn files)
Bashpython TattleTail.py assemblies/ dir_of_ffn/*.ffn -o batch_results --window 15000

Specify custom database path and index prefix
Bashpython TattleTail.py genome.fna -o results \
    --db-fasta database_tailocin.fasta \
    --db-name database_tailocin

Print per-sample reports to terminal (with colors) + export hit details
Bashpython TattleTail.py genome.fna -o results \
    --per-sample-stdout --force-stdout --dump-hits-tsv

Re-summarize previous results (without re-running analysis)
Bashpython TattleTail.py collect \
    --results-root results \
    --out-tsv all_results.tsv \
    --out-json all_results.json \
    --print

For full command-line options:
Bashpython TattleTail.py --help

Database Preparation
TattleTail requires a custom BLAST protein database (database_tailocin.fasta).
Required header format:
text>ID LABEL [optional description]
Example:
text>BAR70105.1 integrase_1 [Pseudomonas aeruginosa]

Build the database (run only once):
Bashmakeblastdb -in database_tailocin.fasta -dbtype prot -out database_tailocin -parse_seqids
Important: The reference database file is not included in this repository due to size and licensing considerations.
Please prepare your own database or contact the author for details.

Output Files
Per-sample outputs (in each sample's result folder)

query_proteins.faa — Predicted or translated proteins
blast_results.txt — Raw BLASTP output
tailocin_report.txt — Human-readable summary report
tailocin_report.json — Structured JSON report
run.log — Execution log
cluster_annotation.tsv — Detailed hit annotations for candidate clusters (if any)

Batch summary (in the root output directory)
allresultsofN samples_YYYYMMDD_HHMM.txt — Plain text combined report
allresultsofN samples_YYYYMMDD_HHMM_colored.txt — Colored ANSI version (view with less -R or cat)













