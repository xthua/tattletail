# TattleTail: Pyocin / Tailocin Prediction Tool (v1.0)

TattleTail is a Python-based bioinformatics pipeline for identifying potential **tailocin / pyocin gene clusters** in bacterial genomes using BLAST-based functional annotation and window-based cluster evaluation.

The tool supports both annotated and unannotated genome inputs and performs automatic gene prediction, BLAST search, clustering, and candidate validation.

---

## Features

- Automatic BLAST database self-check and validation
- Supports mixed input formats (FFN or unannotated FASTA)
- Gene prediction using Prodigal or direct translation from FFN
- BLASTP search with functional label mapping
- Window-based clustering of hits
- Strict four-step candidate evaluation
- Detailed per-cluster annotation output
- Batch processing support
- Plain text + colored ANSI summary reports
- Result re-collection without re-running analysis

---

## Installation

Clone the repository:

```bash
git clone https://github.com/xthua/tattletail.git
cd tattletail
```

Create and activate a conda environment:

```bash
conda create -n tailocin python=3.9
conda activate tailocin
```

Install dependencies:

```bash
conda install -c bioconda prodigal blast biopython
pip install termcolor
```

(If `termcolor` is already installed via conda, the pip step can be skipped.)

---

### Create environment from file

Instead of installing packages individually, you can create the full environment with:

```bash
conda env create -f environment.yml
conda activate tailocin_env
```

---

## Requirements

### Python packages

- Python ≥ 3.9
- Biopython ≥ 1.85
- termcolor ≥ 2.1.0

### External software (must be in PATH)

- Prodigal ≥ 2.6.3
- BLAST+ ≥ 2.16.0

Required BLAST programs:

- blastp
- makeblastdb
- blastdbcmd

Check installation:

```bash
prodigal -v
blastp -version
makeblastdb -version
blastdbcmd -version
```

---

## Usage

### Single sample (FASTA input)

15 kb window, minimum cluster span ≥ 13.4 kb

```bash
python TattleTail.py genome.fna \
    -o results \
    --window 15000 \
    --min-cluster-span 13400
```

---

### Batch processing (mixed .fna / .ffn)

```bash
python TattleTail.py assemblies/*.fna dir_of_ffn/*.ffn \
    -o batch_results \
    --window 15000
```

---

### Specify custom database

```bash
python TattleTail.py genome.fna \
    -o results \
    --db-fasta database_tailocin.fasta \
    --db-name database_tailocin
```

---

### Print per-sample reports + export hit details

```bash
python TattleTail.py genome.fna \
    -o results \
    --per-sample-stdout \
    --force-stdout \
    --dump-hits-tsv
```

---

### Re-summarize previous results (no re-run)

```bash
python TattleTail.py collect \
    --results-root results \
    --out-tsv all_results.tsv \
    --out-json all_results.json \
    --print
```

---

### Show all options

```bash
python TattleTail.py --help
```

---

## Database preparation

TattleTail requires a custom BLAST protein database.

Required FASTA header format:

```
>ID LABEL [optional description]
```

Example:

```
>BAR70105.1 integrase_1 [Pseudomonas aeruginosa]
```

Build database:

```bash
makeblastdb \
    -in database_tailocin.fasta \
    -dbtype prot \
    -out database_tailocin 
```

---

## Output

### Per-sample output folder

Each sample produces:

```
query_proteins.faa
blast_results.txt
tailocin_report.txt
tailocin_report.json
run.log
cluster_annotation.tsv
```

Description:

| File | Description |
|------|------------|
| query_proteins.faa | predicted / translated proteins |
| blast_results.txt | raw BLASTP output |
| tailocin_report.txt | readable report |
| tailocin_report.json | structured report |
| run.log | execution log |
| cluster_annotation.tsv | detailed cluster hits |

---

### Batch summary

In root output directory:

```
allresultsofN samples_YYYYMMDD_HHMM.txt
allresultsofN samples_YYYYMMDD_HHMM_colored.txt
```

Colored file can be viewed with:

```bash
less -R file.txt
```

---

## Citation

If you use **TattleTail** in your research, please cite our preprint:

**TattleTail: A Pyocin / Tailocin Prediction Tool**  
Rayhaan G. Pais, Weilian Chen, Sebastian Leptihn, Xiaoting Hua and Belinda Loh(2026)  
bioRxiv  
DOI: [https://doi.org/10.64898/2026.03.25.712926](https://doi.org/10.64898/2026.03.25.712926)


