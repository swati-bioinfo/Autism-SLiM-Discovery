"""
Phase 7: Fisher's exact test enrichment for MEME-discovered motifs.

Reads MEME XML, extracts motif PWMs, scans autism and control IDR sets,
builds 2x2 contingency tables, runs one-sided Fisher's exact test,
applies Benjamini-Hochberg FDR correction, calculates effect sizes.
"""

import xml.etree.ElementTree as ET
import numpy as np
from scipy.stats import fisher_exact
import csv
import os

MEME_XML = os.path.join("03-analysis", "meme.xml")
AUTISM_IDRS = os.path.join("00-raw", "autism-idrs-clean.fasta")
CONTROL_IDRS = os.path.join("01-control-data", "control-idrs-clean.fasta")
OUTPUT = os.path.join("03-analysis", "enrichment-results.csv")

LETTERS_ORDER = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
letter_to_idx = {l: i for i, l in enumerate(LETTERS_ORDER)}

# --- Parse MEME XML ---
tree = ET.parse(MEME_XML)
root = tree.getroot()

bg = {}
for bf in root.findall(".//background_frequencies"):
    if bf.get("source") == "--sequences--":
        for val in bf.findall(".//value"):
            bg[val.get("letter_id")] = float(val.text)
        break

bg_arr = np.array([bg[l] for l in LETTERS_ORDER])
eps = 1e-10
bg_arr_safe = np.maximum(bg_arr, eps)

motifs = []
for motif_elem in root.findall(".//motif"):
    m = {
        "id": motif_elem.get("id"),
        "name": motif_elem.get("name"),
        "width": int(motif_elem.get("width")),
        "sites": int(motif_elem.get("sites")),
        "e_value": float(motif_elem.get("e_value")),
        "bayes_threshold": float(motif_elem.get("bayes_threshold")),
    }
    ppm_rows = []
    prob_elem = motif_elem.find("probabilities")
    for array in prob_elem.findall("alphabet_matrix/alphabet_array"):
        probs = {}
        for val in array.findall("value"):
            probs[val.get("letter_id")] = float(val.text)
        ppm_rows.append([probs[l] for l in LETTERS_ORDER])
    ppm = np.array(ppm_rows)
    ppm_safe = np.maximum(ppm, eps)
    m["pwm"] = np.log2(ppm_safe / bg_arr_safe[np.newaxis, :])
    motifs.append(m)

print(f"Parsed {len(motifs)} motifs from MEME XML")

# --- Read FASTA ---
def read_fasta(path):
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            seqs[current_id] = "".join(current_seq)
    return seqs

autism_seqs = read_fasta(AUTISM_IDRS)
control_seqs = read_fasta(CONTROL_IDRS)
na, nc = len(autism_seqs), len(control_seqs)
print(f"Autism IDRs: {na} | Control IDRs: {nc}")
print()

# --- Scan sequences with PWM ---
def scan_sequences(seq_dict, pwm, bayes_threshold):
    width = pwm.shape[0]
    count = 0
    for seq in seq_dict.values():
        if len(seq) < width:
            continue
        best = -np.inf
        for i in range(len(seq) - width + 1):
            idxs = [letter_to_idx.get(c) for c in seq[i:i+width]]
            if any(idx is None for idx in idxs):
                continue
            score = pwm[np.arange(width), idxs].sum()
            if score > best:
                best = score
        if best >= bayes_threshold:
            count += 1
    return count

results = []
for m in motifs:
    a = scan_sequences(autism_seqs, m["pwm"], m["bayes_threshold"])
    b = scan_sequences(control_seqs, m["pwm"], m["bayes_threshold"])
    c = na - a
    d = nc - b

    odds_ratio, p_val = fisher_exact([[a, b], [c, d]], alternative="greater")
    enrich_fold = (a / na) / (b / nc) if b > 0 else float("inf")

    results.append({
        "motif_id": m["id"],
        "motif_name": m["name"],
        "width": m["width"],
        "e_value": m["e_value"],
        "bayes_threshold": round(m["bayes_threshold"], 2),
        "autism_with": a, "autism_without": c,
        "control_with": b, "control_without": d,
        "odds_ratio": round(odds_ratio, 4),
        "p_value": p_val,
        "enrichment_fold": round(enrich_fold, 2),
    })
    print(f"{m['id']:>8} ({m['name']}): "
          f"Autism {a}/{na} ({a/na*100:.1f}%) | "
          f"Control {b}/{nc} ({b/nc*100:.1f}%) | "
          f"OR={odds_ratio:.4f} p={p_val:.4e} fold={enrich_fold:.2f}")

# --- Benjamini-Hochberg FDR ---
p_vals = np.array([r["p_value"] for r in results])
n = len(p_vals)
sorted_idx = np.argsort(p_vals)
q_temp = p_vals[sorted_idx] * n / np.arange(1, n + 1)
q_sorted = np.minimum.accumulate(q_temp[::-1])[::-1]
q_values = np.zeros(n)
for i, idx in enumerate(sorted_idx):
    q_values[idx] = q_sorted[i]

for i, r in enumerate(results):
    r["q_value"] = round(q_values[i], 6)
    r["significant_q05"] = q_values[i] < 0.05

print(f"\n--- Benjamini-Hochberg correction (q < 0.05) ---")
for r in results:
    flag = " ***" if r["significant_q05"] else ""
    print(f"  {r['motif_id']} (E={r['e_value']:.2e}): p={r['p_value']:.4e}  q={r['q_value']:.6f}{flag}")

# --- Write CSV ---
fieldnames = [
    "motif_id", "motif_name", "width", "e_value", "bayes_threshold",
    "autism_with", "autism_without", "control_with", "control_without",
    "odds_ratio", "p_value", "q_value", "enrichment_fold",
    "significant_q05",
]
with open(OUTPUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    w.writeheader()
    w.writerows(results)

print(f"\nResults saved to {OUTPUT}")
