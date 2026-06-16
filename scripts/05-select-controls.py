"""Select length-matched control proteins from human Swiss-Prot (non-SFARI)."""
import random
from Bio import SeqIO

random.seed(42)

sfari_records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))
sfari_accs = {r.id.split("|")[1] if "|" in r.id else r.id for r in sfari_records}

human_records = list(SeqIO.parse("00-raw/uniprot_sprot_human.fasta", "fasta"))
candidates = []
for r in human_records:
    acc = r.id.split("|")[1] if "|" in r.id else r.id
    if acc not in sfari_accs:
        candidates.append((len(r.seq), r))

candidates.sort(key=lambda x: x[0])

def find_match(target_len, tol=0.10):
    lo = int(target_len * (1 - tol))
    hi = int(target_len * (1 + tol))
    pool = [r for l, r in candidates if lo <= l <= hi]
    while not pool and tol < 0.5:
        tol += 0.05
        lo = int(target_len * (1 - tol))
        hi = int(target_len * (1 + tol))
        pool = [r for l, r in candidates if lo <= l <= hi]
    if not pool:
        return None
    return random.choice(pool)

matched = []
for sfari_rec in sfari_records:
    target_len = len(sfari_rec.seq)
    match = find_match(target_len)
    if match:
        matched.append(match)

SeqIO.write(matched, "01-control-data/control-proteins.fasta", "fasta")

print(f"SFARI proteins: {len(sfari_records)}")
print(f"Controls selected: {len(matched)}")

sfari_lens = sorted(len(r.seq) for r in sfari_records)
ctrl_lens = sorted(len(r.seq) for r in matched)
print(f"SFARI lengths: min={min(sfari_lens)}, max={max(sfari_lens)}, median={sfari_lens[len(sfari_lens)//2]}")
print(f"Control lengths: min={min(ctrl_lens)}, max={max(ctrl_lens)}, median={ctrl_lens[len(ctrl_lens)//2]}")
