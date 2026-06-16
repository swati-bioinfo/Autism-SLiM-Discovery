# LR-0003: User Exercise — Autism Protein Motif Discovery Project



**Topic:** Autism-linked protein motif discovery

**Comparison:** Conserved short linear motifs in disordered regions of autism risk proteins vs motifs in a matched set of random human proteins.

**Data:** SFARI Gene database + UniProt protein sequences + MEME Suite (web-based motif discovery tool).

**Refined Question:** "Do intrinsically disordered regions of proteins encoded by high-confidence autism risk genes contain conserved short linear motifs that are significantly enriched compared to the rest of the human proteome?"

**FINER Check:** All Pass. Feasible (public data, no code needed), Interesting (shared mechanisms), Novel (underexplored), Ethical (public data), Relevant (biomarkers/drug targets).

**Hypothesis:** "The disordered regions of high-confidence autism risk proteins contain at least one short linear motif (6–15 amino acids) that is significantly enriched (p < 0.05 after multiple testing correction) compared to a length-matched set of random human proteins, as detected by MEME motif discovery."


**Design Brief:**

**Q1 — Controls:**
- **Background:** A matched set of non-autism human proteins (length-matched, from the same proteome) as the baseline for enrichment comparison.
- **Negative control:** Randomly shuffled sequences preserving amino acid composition. If MEME finds "motifs" in shuffled sequences, those are false positives from composition bias.
- **Positive control:** A known short linear motif (e.g., SH3-binding PxxP motif) in a protein known to contain it. If the pipeline misses it, something is wrong.

**Q2 — Replicates (Robustness):**
- **Algorithmic replicates:** Run MEME multiple times with different random seeds. Do the same motifs appear consistently?
- **Method replicates:** Use a second tool (e.g., DREME or SLiMProb) and check if top motifs overlap.
- **Bootstrap replicates:** Randomly sample 80% of autism proteins 100 times. Which motifs appear in ≥ 90% of subsets? Those are robust signals.

**Q3 — Confounders:**
- **Protein length bias:** Autism risk proteins may be longer → more IDR sequence → more chance of finding motifs. *Solution:* Length-match the control set.
- **Sequence composition bias:** IDRs are enriched in P, E, S, Q. Any motif found might reflect composition, not biology. *Solution:* Shuffled-sequence negative control.
- **Annotation bias:** Autism risk proteins may be more studied → more known motifs already annotated. *Solution:* Use de novo motif discovery (MEME) that doesn't rely on existing annotations.

**Q4 — Data Management Plan:**
- **Directory structure:**
  ```
  autism-motif-project/
  ├── 00-raw/            # SFARI gene list, UniProt sequences (read-only)
  ├── 01-control-data/   # Length-matched control proteins, shuffled decoys
  ├── 02-motif-discovery/ # MEME input/output, DREME results
  ├── 03-analysis/       # Enrichment tests, bootstrap results, figures
  └── 04-docs/           # README, parameter log, design brief
  ```
- **Parameter tracking:** Maintain a plain-text `04-docs/parameters.log` where every run records: tool name + version, all non-default parameters, random seed, input file, date. Example:
  ```
  2026-06-16  MEME v5.5.2  --dna --mod zoops --nmotifs 5 --minw 6 --maxw 15 --evt 0.05  seed=42  input=autism_idrs.fa
  ```
  One entry per run. Satisfies Sandve Rules 1 (track provenance) and 3 (archive external program versions).

---

## Data Acquisition Exercise 

**Date of lookup:** 2026-06-16

**Q1 — SFARI Gene:**
- The high-confidence gene list is the **Category 1** list within the **Human Gene** scoring module.
- **Category 1 (High Confidence): 245 genes.**
- Other categories: S (Syndromic) — 218 genes, Category 2 — 710 genes, Category 3 — 228 genes.
- Latest release: **2026 Q1**. Database updated: **May 1, 2026**.
- Total scored genes: 1,183 (94 uncategorized).

**Q2 — UniProt (SHANK3):**
- Reviewed (Swiss-Prot) entry accession: **Q9BYB0**.
- Protein name: SH3 and multiple ankyrin repeat domains protein 3 (SHAN3_HUMAN).
- Length: **1,806 amino acids**.
- Gene names: SHANK3, KIAA1650, PROSAP2, PSAP2.
- UniProt release: **2026_02**.
- *Note: Entry page would need to be checked for specific disordered region / compositional bias annotations under "Region" or "Compositional bias" sections.*

**Q3 — MEME Suite:**
- Version at bottom: **MEME Suite 5.5.9** (previous version 5.5.8 also available via link).
- MEME stands for **Multiple Em for Motif Elicitation**.
- Input format required: **FASTA** (documented at `doc/fasta-format.html`).
- Site posted/updated: **4/6/2026**.

**Q4 — Metadata for reproducibility:**
For each downloaded file, the following metadata should be recorded so another researcher can exactly reproduce the acquisition:

| Field | Example |
|---|---|
| Source URL | `https://gene.sfari.org/database/human-gene/` |
| Download date | `2026-06-16` |
| Database version | `2026 Q1` (SFARI) / `2026_02` (UniProt) |
| File name | `sfari-category-1-genes.csv` |
| Search query | `SHANK3 human` filtered to Swiss-Prot, human |
| Tool used | Manual download via browser / `wget` / API script |
| File format | CSV (SFARI), FASTA (UniProt) |
| Checksum (if available) | `sha256sum filename` |

**Key insight:** Recording metadata is the data-acquisition equivalent of tracking parameters in an experiment — it satisfies Sandve Rule 1 (track provenance) and ensures another researcher can audit exactly what was downloaded.

---

## QC Checklist Exercise 

**Q1 — SFARI Version:**
- Active release at time of lookup: **2026 Q1**, database updated **May 1, 2026**.
- To check for changes: revisit `gene.sfari.org` and look for the "Latest Release" badge in the top navigation bar and the "Database updated on" date on the Gene Scoring page.

**Q2 — Sequence Sources (Swiss-Prot vs TrEMBL):**
- Use **Swiss-Prot (reviewed) only** for maximum confidence — these entries are expert-curated and manually verified.
- If some Category 1 genes lack Swiss-Prot entries: either (a) exclude them and document the gap as a limitation, or (b) use the best TrEMBL entry but flag it explicitly in the QC report.
- Most SFARI Category 1 genes should have Swiss-Prot entries; a preliminary count of missing entries is the first QC step.

**Q3 — IDR Prediction Tools:**
Two widely used tools to fill gaps where UniProt does not annotate disordered regions:
1. **IUPred3** (`iupred3.elte.hu`) — the most widely used IDR predictor, estimates disorder tendency from amino acid pairwise energies.
2. **ESPritz** — a fast alternative using bidirectional recursive neural networks.
Both accept FASTA input and output per-residue disorder scores.

**Q4 — Control Matching QC Checks:**
Three checks to confirm the control set is valid:
1. **Length distribution:** Plot histograms of autism vs control IDR lengths — they should overlap closely.
2. **Amino acid composition:** Compare per-residue frequencies between sets. Major deviations in IDR-enriched residues (P, E, S, Q) signal a biased control.
3. **No overlap:** Confirm no control gene appears in the SFARI autism risk list.

**Q5 — QC Documentation Plan:**
"I will maintain a single `04-docs/qc-report.md` file with four sections — one per dataset (SFARI gene list, autism UniProt sequences, control sequences, shuffled negative controls) — each recording source URL, download date, version, file format, sequence count, length statistics (min, max, mean, median), character validation results, and any filtered records. The report will also include a composition comparison table and a log of every filtering decision with justification. This satisfies Sandve Rules 1 (provenance) and 5 (record intermediate results)."

---

## Analysis & Publication Plan Exercise 

**Q1 — MEME Parameters:**

**Motif width range: 6–15 amino acids.**
- *Minimum 6:* SLiMs shorter than 6 residues (e.g., 3–5 aa) are too short to be statistically discriminative — they occur frequently by chance in any sequence, leading to excessive false positives. A 6-residue motif provides enough specificity to distinguish signal from background.
- *Maximum 15:* While some linear motifs can be longer, the classic SLiM definition caps at ~15 residues. Longer motifs approach domain-level length and would require different discovery strategies (e.g., HMM profiles). IDR-mediated interactions are typically short and modular.

**Model type: zoops (Zero or One Occurrence Per Sequence).**
- *Why zoops:* Each autism IDR sequence likely contains at most one instance of a given SLiM. Using `anops` (Any Number of Occurrences Per Sequence) would inflate false positives by allowing MEME to fit multiple weak instances. Using `oops` (Exactly One Per Sequence) is too restrictive — a motif may genuinely be absent from some sequences.
- *zoops is MEME's default model for protein motif discovery for good reason:* it matches the biological expectation that functional SLiMs are present in a subset of sequences rather than universally.

**Additional parameters:**
- Number of motifs to find: `--nmotifs 10` (cast a wide net; the statistical filter in Phase 6 will prune false positives).
- E-value threshold: `--evt 0.05` (stringent to avoid MEME reporting hundreds of weak motifs).

**Q2 — Statistical Testing Plan for 8 Candidate Motifs:**

For each of the 8 motifs discovered by MEME, I will construct a 2x2 contingency table:

| | Autism IDRs | Control IDRs |
|---|---|---|
| Motif present | a | b |
| Motif absent | c | d |

- **Test:** Fisher's exact test (one-sided, testing enrichment in autism > control). Fisher's is preferred over chi-square when expected cell counts are low (< 5), which is common for rare motifs.
- **Multiple testing correction:** Benjamini-Hochberg False Discovery Rate (FDR) applied to all 8 p-values. This is appropriate because (a) 8 tests is a modest multiple testing burden, and (b) FDR controls the expected proportion of false positives rather than the family-wise error rate, giving more power to detect real motifs.
- **Significance threshold:** q < 0.05 (adjusted p-value).
- **Effect size reporting:** For each significant motif, report the odds ratio (a*d / b*c) and the enrichment fold (prevalence in autism / prevalence in control), not just the p-value.

**Q3 — Validation Strategies for a Motif Passing FDR:**

Two strategies I would apply before claiming biological significance:

1. **Bootstrap resampling (stability test):** Randomly sample 80% of the autism IDR set (with replacement) 1,000 times. For each resample, re-run MEME and check whether the motif is rediscovered. If the motif appears in ≥ 90% of resamples, it is robust to sampling variation. This directly implements the bootstrap replication plan from Lesson 3 (Experimental Design — Q2).

2. **TOMTOM search against ELM database (biological contextualization):** Submit the motif's position-specific scoring matrix (PSSM) to TOMTOM, comparing against the ELM (Eukaryotic Linear Motif) database. A match to a known functional SLiM (e.g., SH3-binding PxxP, 14-3-3 binding mode, nuclear localization signal) provides independent evidence that the discovered pattern is biologically meaningful. A motif with no ELM match is not necessarily spurious — it could be novel — but requires stronger validation before publication.

**Q4 — Strongest Claim the Evidence Could Support:**

"We found evidence that an [X-residue motif with consensus sequence Y] is significantly enriched in the intrinsically disordered regions of proteins encoded by high-confidence autism risk genes (SFARI Category 1) compared to length-matched controls, and this motif matches known [biological function, e.g., SH3-binding / phosphorylation / 14-3-3 binding] motifs in the ELM database."

*This statement is bounded: it specifies the exact dataset (SFARI Category 1), the exact comparison (length-matched controls), the statistical criterion (significant enrichment), and the biological annotation (ELM match). It does not claim causality, mechanism, or therapeutic relevance — those would require additional experiments.*

**Q5 — Three Reproducibility Actions:**

| # | Action | Sandve Rule |
|---|---|---|
| 1 | **Archive all MEME execution parameters** in `04-docs/parameters.log`: MEME version (5.5.9), motif width range, model type, random seed, E-value threshold, and input file checksum. This ensures the exact discovery run can be reproduced. | Rule 3 (Archive external program versions) + Rule 6 (Note random seeds) |
| 2 | **Version-control all analysis code** (scripts for IDR extraction, statistical tests, figure generation) in a Git repository committed to GitHub. Use meaningful commit messages that reference specific analysis steps. Include a `README.md` that documents exactly how to reproduce each figure. | Rule 4 (Version control all custom scripts) |
| 3 | **Deposit raw data** (SFARI gene list, UniProt sequences, control sequences) in Zenodo with a DOI. Include the QC report, the complete MEME output (all motifs found, not just significant ones), and a reproducibility statement in the repository README. | Rule 10 (Provide public access to scripts, runs, and results) |


