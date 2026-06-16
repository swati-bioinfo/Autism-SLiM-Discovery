"""Fetch remaining control IUPred3 predictions with threading (faster)."""
import time, requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO

THRESHOLD = 0.5
MIN_LEN = 15

records = list(SeqIO.parse("01-control-data/control-proteins.fasta", "fasta"))

processed_accs = set()
with open("01-control-data/control-idrs.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) >= 2:
                processed_accs.add(parts[1])

pending = [r for r in records if (r.id.split("|")[1] if "|" in r.id else r.id) not in processed_accs]
print(f"Processing {len(pending)} remaining proteins...")

session = requests.Session()

def fetch(rec):
    acc = rec.id.split("|")[1] if "|" in rec.id else rec.id
    seq_str = str(rec.seq)
    try:
        resp = session.get(f"https://iupred3.elte.hu/iupred3/{acc}", timeout=30)
        if resp.status_code != 200:
            return (acc, None, None)
        scores = []
        for line in resp.text.splitlines():
            if line.startswith("#") or not line.strip(): continue
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                scores.append(float(parts[2]))
        if len(scores) != len(seq_str):
            return (acc, None, f"length mismatch {len(scores)} vs {len(seq_str)}")
        return (acc, (rec, scores), None)
    except Exception as e:
        return (acc, None, str(e))

new_idrs = 0
errors = 0
with open("01-control-data/control-idrs.fasta", "a") as out, ThreadPoolExecutor(max_workers=5) as pool:
    futures = [pool.submit(fetch, rec) for rec in pending]
    for i, fut in enumerate(as_completed(futures)):
        acc, result, err = fut.result()
        if err or result is None:
            errors += 1
            print(f"  ERROR [{acc}]: {err}")
            continue
        rec, scores = result
        seq_str = str(rec.seq)
        in_idr = False
        start = 0
        for j in range(len(scores)):
            if scores[j] > THRESHOLD and not in_idr:
                start = j
                in_idr = True
            elif scores[j] <= THRESHOLD and in_idr:
                if j - start >= MIN_LEN:
                    out.write(f">{rec.id}_IDR_{start+1}-{j}\n{seq_str[start:j]}\n")
                    new_idrs += 1
                in_idr = False
        if in_idr and len(scores) - start >= MIN_LEN:
            out.write(f">{rec.id}_IDR_{start+1}-{len(scores)}\n{seq_str[start:]}\n")
            new_idrs += 1
        time.sleep(0.1)

print(f"\nDone. Added {new_idrs} new IDRs, {errors} errors.")
