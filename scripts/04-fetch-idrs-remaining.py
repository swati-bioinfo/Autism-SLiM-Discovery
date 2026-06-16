"""Process remaining proteins not yet in autism-idrs.fasta."""
import time, requests
from Bio import SeqIO

THRESHOLD = 0.5
MIN_LEN = 15

all_recs = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))

processed = set()
with open("00-raw/autism-idrs.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) >= 2:
                processed.add(parts[1])

remaining = [r for r in all_recs if (r.id.split("|")[1] if "|" in r.id else r.id) not in processed]
print(f"Remaining to process: {len(remaining)}")

idr_count = 0
with open("00-raw/autism-idrs.fasta", "a") as out:
    for i, rec in enumerate(remaining):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_str = str(rec.seq)

        try:
            resp = requests.get(f"https://iupred3.elte.hu/iupred3/{acc}", timeout=30)
        except:
            print(f"  [{acc}] failed")
            time.sleep(0.5)
            continue

        if resp.status_code != 200 or not resp.text.strip():
            print(f"  [{acc}] bad response")
            time.sleep(0.5)
            continue

        scores = []
        for line in resp.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                scores.append(float(parts[2]))

        if len(scores) != len(seq_str):
            print(f"  [{acc}] length mismatch")
            time.sleep(0.5)
            continue

        in_idr = False
        start = 0
        for j in range(len(scores)):
            if scores[j] > THRESHOLD and not in_idr:
                start = j
                in_idr = True
            elif scores[j] <= THRESHOLD and in_idr:
                length = j - start
                if length >= MIN_LEN:
                    out.write(f">{rec.id}_IDR_{start+1}-{j}\n{seq_str[start:j]}\n")
                    idr_count += 1
                in_idr = False
        if in_idr:
            length = len(scores) - start
            if length >= MIN_LEN:
                out.write(f">{rec.id}_IDR_{start+1}-{len(scores)}\n{seq_str[start:]}\n")
                idr_count += 1

        time.sleep(0.5)
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(remaining)} — {idr_count} new IDRs")

print(f"Done. Added {idr_count} new IDR regions.")
