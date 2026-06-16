"""Validate fetched sequences: count, lengths, characters."""
from Bio import SeqIO

records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))
print(f"Total sequences: {len(records)}")

lengths = [len(r.seq) for r in records]
print(f"Length range: {min(lengths)} - {max(lengths)} AA")
print(f"Mean length: {sum(lengths)/len(lengths):.0f} AA")

valid = set("ACDEFGHIKLMNPQRSTVWY")
bad = []
for r in records:
    for c in str(r.seq).upper():
        if c not in valid:
            bad.append((r.id, c))
if bad:
    print(f"Found {len(bad)} invalid characters (e.g., B, Z, X)")
else:
    print("All sequences use standard 20 amino acids -- clean!")

with open("04-docs/qc-summary.txt", "w") as f:
    f.write(f"Sequences: {len(records)}\n")
    f.write(f"Min length: {min(lengths)}\n")
    f.write(f"Max length: {max(lengths)}\n")
    f.write(f"Mean length: {sum(lengths)/len(lengths):.0f}\n")
    f.write(f"Invalid chars: {len(bad)}\n")
