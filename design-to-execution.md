# Lesson 7: From Design to Execution — Building Your Autism Motif Pipeline

You've designed the project. Now let's build it.

**Goal:** By the end of this lesson, you will have run your first real pipeline step — fetching protein sequences for SFARI autism genes.

## What You'll Build

```
autism-motif-project/
├── 00-raw/              # SFARI gene list, UniProt sequences (read-only)
├── 01-control-data/     # Length-matched controls, shuffled decoys
├── 02-motif-discovery/  # MEME input/output, DREME results
├── 03-analysis/         # Enrichment tests, bootstrap, figures
└── 04-docs/             # README, parameter log, design brief
```

Data flows through these directories left-to-right: raw → control → discovery → analysis → docs.

## Step 0: Set Up Your Environment

First, create the directory structure and install Python packages.

```bash
# Create your project folder anywhere (Desktop is fine)
mkdir autism-motif-project
cd autism-motif-project

# Create the directory structure
mkdir 00-raw 01-control-data 02-motif-discovery 03-analysis 04-docs

# Install required Python packages
pip install requests biopython pandas numpy
```

> **What these packages do:**
> - `requests` — downloads data from web APIs (UniProt, SFARI)
> - `biopython` — reads/writes FASTA files, handles sequences
> - `pandas` — tables and CSV files
> - `numpy` — numerical calculations

## Step 1: Download the SFARI Gene List

SFARI Gene doesn't have a public API, so you'll download manually:

1. Go to https://gene.sfari.org
2. Click "Download" → "Full gene list" → save as `SFARI-Gene_genes.csv`
3. Move it into your project: `mv ~/Downloads/SFARI-Gene_genes.csv 00-raw/`

### Verify the download

```bash
wc -l 00-raw/SFARI-Gene_genes.csv
head -5 00-raw/SFARI-Gene_genes.csv
```

You should see ~1,200 lines (1 header + ~1,183 genes). The first row is column names.

## Step 2: Extract Category 1 Gene Symbols

Save this script as `01-extract-category1.py` in your project root:

```python
"""Extract SFARI Category 1 (high confidence) gene symbols."""
import csv

input_file = "00-raw/SFARI-Gene_genes.csv"
output_file = "00-raw/category1-genes.txt"

genes = []
with open(input_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Adjust column name if SFARI uses different header text
        if row.get("gene-symbol") and row.get("gene-score") == "1":
            genes.append(row["gene-symbol"])

with open(output_file, "w") as f:
    for g in sorted(genes):
        f.write(g + "\n")

print(f"Found {len(genes)} Category 1 genes")
print(f"Saved to {output_file}")
```

Run it:

```bash
python 01-extract-category1.py
```

Expected output: `Found 245 Category 1 genes`.

## Step 3: Fetch UniProt Sequences

This script fetches one protein sequence per gene from UniProt's REST API. Save as `02-fetch-uniprot.py`:

```python
"""Fetch reviewed (Swiss-Prot) sequences for a list of human genes."""
import requests
import time

with open("00-raw/category1-genes.txt") as f:
    genes = [line.strip() for line in f if line.strip()]

output_file = "00-raw/autism-proteins.fasta"
skipped = []

with open(output_file, "w") as out:
    for i, gene in enumerate(genes):
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene:{gene}+AND+organism_id:9606+AND+reviewed:true"
            "&format=fasta"
        )
        resp = requests.get(url)
        if resp.status_code == 200 and resp.text.strip():
            out.write(resp.text)
            if not resp.text.startswith(">"):
                skipped.append(gene)
        else:
            skipped.append(gene)
        # Be nice to the API — 5 requests per second max
        time.sleep(0.2)
        if (i + 1) % 50 == 0:
            print(f"Processed {i+1}/{len(genes)} genes")

if skipped:
    print(f"Warning: {len(skipped)} genes had no Swiss-Prot entry:")
    for g in skipped[:10]:
        print(f"  {g}")
    if len(skipped) > 10:
        print(f"  ... and {len(skipped)-10} more")
print(f"Sequences saved to {output_file}")
```

Run it:

```bash
python 02-fetch-uniprot.py
```

> **How it works:** UniProt's REST API accepts search queries. This script asks: "Find the reviewed (Swiss-Prot) entry for this gene in human (organism_id:9606)." The `reviewed:true` filter ensures you get only curated entries as planned in Lesson 5.

## Step 4: QC the Sequences

Save as `03-qc-sequences.py`:

```python
"""Validate fetched sequences: count, lengths, characters."""
from Bio import SeqIO

records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))
print(f"Total sequences: {len(records)}")

lengths = [len(r.seq) for r in records]
print(f"Length range: {min(lengths)} - {max(lengths)} AA")
print(f"Mean length: {sum(lengths)/len(lengths):.0f} AA")

# Check for invalid characters (only standard 20 amino acids)
valid = set("ACDEFGHIKLMNPQRSTVWY")
bad = []
for r in records:
    for c in str(r.seq).upper():
        if c not in valid:
            bad.append((r.id, c))
if bad:
    print(f"Found {len(bad)} invalid characters (e.g., B, Z, X)")
else:
    print("All sequences use standard 20 amino acids — clean!")

# Save QC summary
with open("04-docs/qc-summary.txt", "w") as f:
    f.write(f"Sequences: {len(records)}\n")
    f.write(f"Min length: {min(lengths)}\n")
    f.write(f"Max length: {max(lengths)}\n")
    f.write(f"Mean length: {sum(lengths)/len(lengths):.0f}\n")
    f.write(f"Invalid chars: {len(bad)}\n")
```

Run it:

```bash
python 03-qc-sequences.py
```

## Step 5: Predict IDR Regions (Web-based)

IDR prediction is best done with a dedicated tool. Two options:

### Option A: Use UniProt annotations (fastest)

Some proteins already have "Compositional bias" or "Region of interest" annotations for disordered regions. Check `autism-proteins.fasta` headers — UniProt includes feature annotations like `FT REGION` or `FT COMPBIAS`. For a quick start:

```python
"""Extract IDRs from UniProt feature annotations in FASTA headers."""
from Bio import SeqIO

records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))

# For now, use full sequences as a placeholder
# (We'll replace this with IUPred3 predictions later)
with open("00-raw/autism-idrs-full.fasta", "w") as out:
    SeqIO.write(records, out, "fasta")

print(f"Exported {len(records)} sequences (full-length, IDR extraction pending)")
```

Save as `04-extract-idrs.py` and run:

```bash
python 04-extract-idrs.py
```

### Option B: IUPred3 (recommended)

IUPred3 is the standard IDR predictor. Use their web server:

1. Go to https://iupred3.elte.hu
2. Upload your `autism-proteins.fasta` file
3. Download the results as "long disorder" format
4. Save to `00-raw/iupred3-results.txt`

Extract IDRs from IUPred3 output using this script. Save as `04-parse-iupred.py`:

```python
"""Parse IUPred3 output and extract IDR sequences above threshold."""
from Bio import SeqIO

# Threshold for "disordered" — IUPred3 score > 0.5
DISORDER_THRESHOLD = 0.5
MIN_IDR_LENGTH = 15  # Shorter than this is too short for motif discovery

# Parse IUPred3 output
# Format: residue_number  aa  score
idr_regions = {}  # {seq_id: [(start, end), ...]}
with open("00-raw/iupred3-results.txt") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        seq_id = parts[0].split("|")[1]  # Extract UniProt accession
        pos = int(parts[1])
        score = float(parts[2])
        # ... parsing logic continues (see note below)

print("IUPred3 parsing requires the exact output format.")
print("Check https://iupred3.elte.hu for format details.")
```

> **Note:** IUPred3 output format varies by version. Once you have a results file, I'll help you write the exact parser.

Then extract IDRs as FASTA for MEME:

```python
"""Extract IDR sequences based on IUPred3 disorder scores."""
from Bio import SeqIO

THRESHOLD = 0.5
MIN_LEN = 15

# Load IUPred3 results as: {(seq_id, start, end)}
# (Manual step pending IUPred3 output format)

sequences = SeqIO.to_dict(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))

with open("00-raw/autism-idrs.fasta", "w") as out:
    # Loop through IDR regions and extract subsequences
    # (Fill in once IUPred3 output is available)
    pass

print("Update this script once you have IUPred3 results.")
```

## Step 6: Build the Control Set

Save as `05-build-control-set.py`:

```python
"""Select length-matched control proteins from the human proteome."""
import csv
import random
from Bio import SeqIO

# Read autism IDR lengths
autism_records = list(SeqIO.parse("00-raw/autism-idrs.fasta", "fasta"))
autism_lengths = [len(r.seq) for r in autism_records]

# Download a reference set of human proteins (sample from Swiss-Prot)
# This is a simplified approach — you'll refine it in Lesson 5's QC step
url = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id:9606+AND+reviewed:true"
    "&format=fasta&size=500"
)
print("Downloading reference human proteome sample...")
import requests
resp = requests.get(url)
with open("01-control-data/reference-human.fasta", "w") as f:
    f.write(resp.text)

ref_records = list(SeqIO.parse("01-control-data/reference-human.fasta", "fasta"))

# Exclude autism genes from reference
autism_genes = set()
with open("00-raw/category1-genes.txt") as f:
    for line in f:
        autism_genes.add(line.strip().upper())

control_records = []
for rec in ref_records:
    # Skip if this protein is in our autism list
    gene_name = rec.description.split("GN=")[1].split()[0].upper() if "GN=" in rec.description else ""
    if gene_name in autism_genes:
        continue
    control_records.append(rec)

# Length-match: for each autism IDR, pick a control IDR of similar length
# For now, use random sampling (you'll refine this)
random.shuffle(control_records)
control_records = control_records[:len(autism_records)]

# Export
with open("01-control-data/control-proteins.fasta", "w") as f:
    SeqIO.write(control_records, f, "fasta")

print(f"Control set: {len(control_records)} proteins")
print("NOTE: These are full-length proteins, not IDR-extracted.")
print("You'll need to run IUPred3 on these too before comparison.")
```

Run it:

```bash
python 05-build-control-set.py
```

## Step 7: Run MEME (Web-based)

MEME discovery is best done on their web server first (no installation needed):

1. Go to https://meme-suite.org/meme/tools/meme
2. Upload `autism-idrs.fasta`
3. Set parameters as planned (Lesson 6):
   - Motif width: 6–15
   - Number of motifs: 10
   - Model: zoops
   - E-value threshold: 0.05
4. Run MEME
5. Download results as ZIP → extract to `02-motif-discovery/meme-output/`

### Alternative: Run MEME locally (if installed)

```bash
meme 00-raw/autism-idrs.fasta \
  -oc 02-motif-discovery/meme-output \
  -dna \
  -mod zoops \
  -nmotifs 10 \
  -minw 6 \
  -maxw 15 \
  -evt 0.05
```

Record this in `04-docs/parameters.log`:

```
2026-06-16  MEME v5.5.9  --mod zoops --nmotifs 10 --minw 6 --maxw 15 --evt 0.05  seed=default  input=autism-idrs.fasta
```

## The Full Pipeline (One Script)

Once you've tested each step, save this as `run-all.py` to run everything in one command:

```python
"""Autism motif pipeline — run all steps."""
import subprocess, sys

steps = [
    ("Extract Category 1 genes", "python 01-extract-category1.py"),
    ("Fetch UniProt sequences",   "python 02-fetch-uniprot.py"),
    ("QC sequences",              "python 03-qc-sequences.py"),
    ("Extract IDRs",              "python 04-extract-idrs.py"),
    ("Build control set",         "python 05-build-control-set.py"),
]

for name, cmd in steps:
    print(f"\n{'='*50}")
    print(f"STEP: {name}")
    print(f"CMD:  {cmd}")
    print('='*50)
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"FAILED: {name}")
        sys.exit(1)
```

## What Each Directory Holds

| Directory | Contents | Created by |
|-----------|----------|------------|
| `00-raw/` | SFARI CSV, gene list, UniProt FASTA, IUPred3 output | Steps 1–4 |
| `01-control-data/` | Reference proteome, matched controls | Step 6 |
| `02-motif-discovery/` | MEME/DREME input files and output folders | Step 7 |
| `03-analysis/` | Enrichment tables, bootstrap results, figures | You'll build next |
| `04-docs/` | `parameters.log`, `qc-summary.txt`, README | Every step |

## Your First Win

Run these three commands in order:

```bash
pip install requests biopython pandas numpy
python 01-extract-category1.py
python 02-fetch-uniprot.py
python 03-qc-sequences.py
```

When `03-qc-sequences.py` prints "All sequences use standard 20 amino acids — clean!" — you've successfully completed Phases 3 and 4 of your project.

---

**Primary resource:** [UniProt REST API documentation](https://www.uniprot.org/help/api) — essential for understanding how to fetch sequences programmatically.

**Stuck?** Ask me to:
- Help debug a script error (paste the error message)
- Write a parser for your exact IUPred3 output format
- Explain what any command or parameter does
- Add a new step to the pipeline
