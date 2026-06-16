"""Autism motif pipeline -- run all automated steps."""
import subprocess
import sys

steps = [
    ("Extract Category 1 genes", "python 01-extract-category1.py"),
    ("Fetch UniProt sequences",  "python 02-fetch-uniprot.py"),
    ("QC sequences",             "python 03-qc-sequences.py"),
    ("Extract IDRs",             "python 04-extract-idrs.py"),
    ("Build control set",        "python 05-build-control-set.py"),
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

print("\nAll steps completed successfully!")
