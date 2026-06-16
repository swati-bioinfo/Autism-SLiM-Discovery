"""Extract full-length sequences as placeholder for IDR extraction.

TODO: Replace this with IUPred3-based IDR parsing once you have
the IUPred3 results file from https://iupred3.elte.hu
"""
from Bio import SeqIO

records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))

with open("00-raw/autism-idrs.fasta", "w") as out:
    SeqIO.write(records, out, "fasta")

print(f"Exported {len(records)} sequences (full-length, IDR extraction pending)")
