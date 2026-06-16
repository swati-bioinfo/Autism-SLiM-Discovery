"""Fetch IUPred3 disorder predictions via REST API, extract IDR sequences."""
import time
import requests
from Bio import SeqIO

DISORDER_THRESHOLD = 0.5
MIN_IDR_LENGTH = 15
API_DELAY = 0.3

records = list(SeqIO.parse("00-raw/autism-proteins.fasta", "fasta"))

idr_count = 0
with open("00-raw/autism-idrs.fasta", "w") as out:
    for i, rec in enumerate(records):
        acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
        seq_str = str(rec.seq)
        url = f"https://iupred3.elte.hu/iupred3/{acc}"

        try:
            resp = requests.get(url, timeout=15)
        except requests.RequestException as e:
            print(f"  [{acc}] API error: {e}")
            time.sleep(API_DELAY)
            continue

        if resp.status_code != 200 or not resp.text.strip():
            print(f"  [{acc}] bad response (status {resp.status_code})")
            time.sleep(API_DELAY)
            continue

        scores = []
        for line in resp.text.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                scores.append(float(parts[2]))

        if len(scores) != len(seq_str):
            print(f"  [{acc}] score/seq length mismatch: {len(scores)} vs {len(seq_str)}")
            time.sleep(API_DELAY)
            continue

        in_idr = False
        start = 0
        for j in range(len(scores)):
            if scores[j] > DISORDER_THRESHOLD and not in_idr:
                start = j
                in_idr = True
            elif scores[j] <= DISORDER_THRESHOLD and in_idr:
                length = j - start
                if length >= MIN_IDR_LENGTH:
                    idr_seq = seq_str[start:j]
                    header = f">{rec.id}_IDR_{start+1}-{j}"
                    out.write(f"{header}\n{idr_seq}\n")
                    idr_count += 1
                in_idr = False
        if in_idr:
            length = len(scores) - start
            if length >= MIN_IDR_LENGTH:
                idr_seq = seq_str[start:]
                header = f">{rec.id}_IDR_{start+1}-{len(scores)}"
                out.write(f"{header}\n{idr_seq}\n")
                idr_count += 1

        if (i + 1) % 20 == 0:
            print(f"Processed {i+1}/{len(records)} — {idr_count} IDRs found")

        time.sleep(API_DELAY)

print(f"\nDone. Processed {len(records)} proteins.")
print(f"Extracted {idr_count} IDR regions to 00-raw/autism-idrs.fasta")
