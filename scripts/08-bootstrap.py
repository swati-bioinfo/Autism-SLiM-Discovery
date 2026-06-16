"""
Bootstrap stability test for MEME Motif 3 (MSTTIMETTTTMATT).

Resamples 80% of autism IDR sequences with replacement (1000x),
scores each against Motif 3 PWM from MEME XML, reports stability.
"""
import sys, math, random, csv
from pathlib import Path
from collections import Counter

BASE = Path(__file__).resolve().parent.parent
XML = BASE / "03-analysis" / "meme.xml"
FASTA = BASE / "00-raw" / "autism-idrs-clean.fasta"
OUT = BASE / "03-analysis" / "bootstrap-results.csv"
SEED = 42
N_ITER = 1000
SAMPLE_FRAC = 0.8

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

def parse_meme_xml(path):
    """Extract Motif 3 PWM (probabilities matrix) from MEME XML."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    # Find motif_3 section
    import re
    # Get the probabilities matrix for motif_3
    m3 = re.search(
        r'<motif id="motif_3".*?</motif>', text, re.DOTALL
    )
    if not m3:
        raise ValueError("motif_3 not found in XML")
    motif_block = m3.group()

    # Extract probabilities matrix (second alphabet_matrix block within motif)
    prob_matches = list(re.finditer(
        r'<alphabet_matrix>\s*((?:<value[^>]*/>\s*)+)</alphabet_matrix>',
        motif_block, re.DOTALL
    ))
    # The second alphabet_matrix in the motif is the probabilities
    # Actually in MEME XML, it's: scores -> probabilities -> regex -> contributions
    # Let's find the probabilities section
    prob_section = re.search(
        r'<probabilities>\s*<alphabet_matrix>(.*?)</alphabet_matrix>\s*</probabilities>',
        motif_block, re.DOTALL
    )
    if not prob_section:
        raise ValueError("probabilities matrix not found for motif_3")

    prob_text = prob_section.group(1)
    # Split into position arrays
    pos_blocks = re.findall(
        r'<alphabet_array>(.*?)</alphabet_array>', prob_text, re.DOTALL
    )
    pwm = []
    for block in pos_blocks:
        values = re.findall(r'<value\s+letter_id="(\w)"[^>]*>([\d.eE+-]+)</value>', block)
        if len(values) != 20:
            # Try alternate format without explicit letter_id
            values = re.findall(r'<value[^>]*>([\d.eE+-]+)</value>', block)
            if len(values) == 20:
                values = list(zip(AA_ORDER, values))
            else:
                raise ValueError(f"Expected 20 values, got {len(values)}")
        probs = {aa: float(v) for aa, v in values}
        pwm.append([probs[aa] for aa in AA_ORDER])
    return pwm  # list of 15 positions, each a list of 20 AA probabilities

def parse_fasta(path):
    """Return list of (seq_id, seq) tuples."""
    seqs = []
    with open(path) as f:
        lines = f.readlines()
    cur_id = None
    cur_seq = []
    for line in lines:
        if line.startswith(">"):
            if cur_id:
                seqs.append((cur_id, "".join(cur_seq)))
            cur_id = line[1:].strip().split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_id:
        seqs.append((cur_id, "".join(cur_seq)))
    return seqs

def score_sequence(seq, pwm):
    """Score a sequence against a PWM. Returns max log-odds score."""
    w = len(pwm)
    if len(seq) < w:
        return -float("inf")
    best = -float("inf")
    bg = 1.0 / 20  # uniform background
    for i in range(len(seq) - w + 1):
        score = 0.0
        for j, aa in enumerate(seq[i:i+w]):
            try:
                aa_idx = AA_ORDER.index(aa)
                prob = pwm[j][aa_idx]
                if prob > 0:
                    score += math.log2(prob / bg)
                else:
                    score += -10  # small penalty for impossible AA
            except ValueError:
                score += -5  # non-standard AA
        if score > best:
            best = score
    return best

def main():
    random.seed(SEED)

    print("Parsing MEME XML for Motif 3 PWM...")
    pwm = parse_meme_xml(XML)
    w = len(pwm)
    print(f"  PWM width: {w}")

    print("Parsing autism IDR FASTA...")
    seqs = parse_fasta(FASTA)
    n_total = len(seqs)
    print(f"  {n_total} sequences")

    print("Scoring all sequences against PWM...")
    scores = [score_sequence(seq, pwm) for _, seq in seqs]
    above_zero = sum(1 for s in scores if s > 0)
    print(f"  {above_zero}/{n_total} sequences score > 0")

    # Determine threshold: the minimum score among the top few sequences
    # MEME found 6 sites. Let's find the 6th highest score.
    sorted_scores = sorted(scores, reverse=True)
    threshold = sorted_scores[min(5, len(sorted_scores)-1)] if sorted_scores else 0
    print(f"  Threshold (6th highest score): {threshold:.2f}")

    # Count sequences in full set that exceed threshold
    n_full = sum(1 for s in scores if s >= threshold)
    print(f"  {n_full}/{n_total} sequences meet threshold in full set")

    print(f"Bootstrapping ({N_ITER} iterations, {SAMPLE_FRAC*100:.0f}% sample)...")
    bootstrap_counts = []
    for it in range(N_ITER):
        sample = random.choices(range(n_total), k=int(n_total * SAMPLE_FRAC))
        count = sum(1 for idx in sample if scores[idx] >= threshold)
        bootstrap_counts.append(count)
        if (it + 1) % 200 == 0:
            print(f"  {it+1}/{N_ITER}")

    mean_c = sum(bootstrap_counts) / N_ITER
    sd_c = (sum((c - mean_c)**2 for c in bootstrap_counts) / N_ITER) ** 0.5
    pct_detected = sum(1 for c in bootstrap_counts if c >= 1) / N_ITER * 100
    pct_ge3 = sum(1 for c in bootstrap_counts if c >= 3) / N_ITER * 100

    print(f"\nResults:")
    print(f"  Mean motif sites per bootstrap: {mean_c:.1f} (SD={sd_c:.1f})")
    print(f"  Motif detected (>=1 site) in: {pct_detected:.1f}% of bootstraps")
    print(f"  Motif with >=3 sites in: {pct_ge3:.1f}% of bootstraps")
    print(f"  Full set count: {n_full}")

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["iteration", "sites_found", "n_total", "n_sampled"])
        for it, c in enumerate(bootstrap_counts):
            w.writerow([it+1, c, n_total, int(n_total * SAMPLE_FRAC)])
        w.writerow([])
        w.writerow(["summary_statistic", "value"])
        w.writerow(["threshold_score", round(threshold, 4)])
        w.writerow(["full_set_count", n_full])
        w.writerow(["mean_bootstrap_count", round(mean_c, 2)])
        w.writerow(["sd_bootstrap_count", round(sd_c, 2)])
        w.writerow(["pct_detected_at_least_1", round(pct_detected, 1)])
        w.writerow(["pct_detected_at_least_3", round(pct_ge3, 1)])
        w.writerow(["n_iterations", N_ITER])
        w.writerow(["sample_fraction", SAMPLE_FRAC])
        w.writerow(["seed", SEED])

    print(f"\nResults saved to {OUT}")
    stability = "STABLE" if pct_detected >= 90 else "BORDERLINE" if pct_detected >= 50 else "UNSTABLE"
    print(f"Verdict: {stability} (detected in {pct_detected:.1f}% of bootstraps)")

if __name__ == "__main__":
    main()
