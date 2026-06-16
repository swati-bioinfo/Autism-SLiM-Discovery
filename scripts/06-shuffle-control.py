"""Generate shuffled IDR sequences as negative control for composition bias.
Shuffles each IDR independently preserving exact AA composition per sequence.
Output: shuffled-autism-idrs.fasta, shuffled-control-idrs.fasta for MEME."""
import random
from Bio import SeqIO

random.seed(99)  # Different seed from the control selection (42)

for label, infile, outfile in [
    ("autism", "00-raw/autism-idrs-clean.fasta", "01-control-data/shuffled-autism-idrs.fasta"),
    ("control", "01-control-data/control-idrs-clean.fasta", "01-control-data/shuffled-control-idrs.fasta"),
]:
    records = list(SeqIO.parse(infile, "fasta"))
    shuffled = []
    for rec in records:
        seq = list(str(rec.seq))
        random.shuffle(seq)
        rec.seq = rec.seq.__class__("".join(seq))
        # Mark as shuffled in the header
        rec.id = rec.id + "_SHUFFLED"
        rec.description = ""
        shuffled.append(rec)
    SeqIO.write(shuffled, outfile, "fasta")
    print(f"{label}: {len(records)} sequences shuffled -> {outfile}")

print(f"\nSeed: 99")
print("Purpose: Negative control for sequence composition bias.")
print("If MEME finds motifs in shuffled sequences, they are composition artifacts.")
