"""Fetch IUPred3 disorder predictions for control proteins via REST API."""
import time, requests
from Bio import SeqIO

THRESHOLD = 0.5
MIN_LEN = 15

records = list(SeqIO.parse("01-control-data/control-proteins.fasta", "fasta"))
total = len(records)

def get_iupred3_scores(accession):
    try:
        resp = requests.get(f"https://iupred3.elte.hu/iupred3/{accession}", timeout=30)
        if resp.status_code != 200:
            return None
        scores = []
        for line in resp.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                scores.append(float(parts[2]))
        return scores
    except:
        return None

idr_count = 0
failed = 0
with open("01-control-data/control-idrs.fasta", "w") as out:
    for i, rec in enumerate(records):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_str = str(rec.seq)

        scores = get_iupred3_scores(acc)
        if scores is None or len(scores) != len(seq_str):
            failed += 1
            print(f"  [{acc}] failed or length mismatch")
            time.sleep(0.3)
            continue

        in_idr = False
        start = 0
        for j in range(len(scores)):
            if scores[j] > THRESHOLD and not in_idr:
                start = j
                in_idr = True
            elif scores[j] <= THRESHOLD and in_idr:
                if j - start >= MIN_LEN:
                    out.write(f">{rec.id}_IDR_{start+1}-{j}\n{seq_str[start:j]}\n")
                    idr_count += 1
                in_idr = False
        if in_idr and len(scores) - start >= MIN_LEN:
            out.write(f">{rec.id}_IDR_{start+1}-{len(scores)}\n{seq_str[start:]}\n")
            idr_count += 1

        if (i + 1) % 30 == 0:
            print(f"  {i+1}/{total} — {idr_count} IDRs found")

        time.sleep(0.3)

print(f"\nDone. {total} proteins, {idr_count} control IDRs extracted, {failed} failed")
