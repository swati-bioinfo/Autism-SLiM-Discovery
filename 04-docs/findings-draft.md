# Findings Draft — Autism SLiM Discovery Project

## Claim 1: Poly-Glutamine Tract (Motif 2) — FDR Significant

> We found evidence that an **11-residue poly-glutamine (poly-Q) motif** is
> significantly enriched in the intrinsically disordered regions of proteins
> encoded by high-confidence autism risk genes (SFARI Category 1, n=242)
> compared to length-matched controls (n=251 non-SFARI human proteins).
>
> **Enrichment:** 41/1,502 autism IDRs (2.7%) vs. 7/1,130 control IDRs (0.6%),
> odds ratio = 4.5, enrichment fold = 4.4×, **FDR q = 0.00022** (Benjamini-Hochberg).
>
> **Interpretation:** Poly-glutamine tracts are a well-known compositional bias
> associated with transcriptional regulation and, when expanded beyond a
> threshold, several neurological disorders (Huntington's disease, SCA types).
> In this autism dataset, **44 unique SFARI genes** contained poly-Q tracts in
> their IDRs, many of which are transcription factors and chromatin regulators
> (e.g., MED13L, KMT2C, KMT2D, CHD2, CHD8, SETD5, ASH1L, TBL1XR1, TCF20).
> The enrichment suggests poly-Q tracts in IDRs may be a shared property of
> autism-risk transcriptional regulators, potentially modulating protein-protein
> interaction networks via homotypic Q/N-rich phase separation.

## Claim 2: Poly-Histidine Tract (Motif 1) — FDR Significant

> We found evidence that a **14-residue poly-histidine motif** (N-terminal Glu,
> consensus EHHHHHHHHHHHHH) is significantly enriched in the intrinsically
> disordered regions of high-confidence autism risk proteins compared to
> length-matched controls.
>
> **Enrichment:** 12/1,502 autism IDRs (0.8%) vs. 0/1,130 control IDRs (0.0%),
> odds ratio = ∞, Fisher's exact p = 0.0012, **FDR q = 0.0059** (Benjamini-Hochberg).
>
> **Interpretation:** Poly-histidine tracts are rare in the human proteome and
> are primarily found in zinc-finger transcription factors and DNA-binding
> proteins. In this set, poly-H tracts appeared in **12 SFARI genes**, most
> notably zinc-finger proteins (ZNFs), TATA-box binding proteins (TBP), and
> transcription factors involved in neuronal development (FOXP1, MEF2C). The'
> complete absence in controls suggests poly-H repeats may be a specific
> signature of autism-related transcriptional machinery, though the low
> absolute count (12 hits) warrants cautious interpretation.

## Claim 3: SPOP-Binding Degron-like Motif (Motif 3) — Bootstrapped Stable

> We found evidence that a **15-residue motif** (consensus
> M[SA]T[TS][IV]METTTT[ML]AT[TS], best ELM match:
> **DEG_SPOP_SBC_1**, SPOP-binding degron, p = 0.005) is enriched in the
> intrinsically disordered regions of high-confidence autism risk proteins
> compared to controls, though the enrichment does not survive FDR correction.
>
> **Enrichment:** 6/1,502 autism IDRs (0.4%) vs. 0/1,130 control IDRs (0.0%),
> odds ratio = ∞, p = 0.034, **FDR q = 0.086 (not significant)**.
>
> **Bootstrap stability (1000×, 80% resample):** Mean sites = 4.8 (SD = 2.2),
> 98.6% of bootstraps detected ≥1 site, 85.0% detected ≥3 sites.
>
> **Interpretation:** Despite falling short of FDR significance, this motif is
> the strongest structured SLiM candidate in our discovery set. The SPOP-binding
> degron match is notable because SPOP (Speckle-type POZ protein) is an E3
> ubiquitin ligase adaptor that targets substrates for proteasomal degradation,
> and SPOP mutations are linked to autism. The six proteins containing this
> motif are strong candidates for SPOP-mediated regulation in an autism context.
> The motif is robustly detected under resampling, indicating it is not driven
> by outlier sequences.

## Summary Table

| Motif | Consensus | Width | E-value | OR | q-value | FDR sig? | Bootstrap stable? | ELM hit |
|-------|-----------|-------|---------|-----|---------|----------|-----------------|---------|
| 1 (poly-H) | EHHHHHHHHHHHHH | 14 | 3.1e-31 | ∞ | 0.0059 | **Yes** | — | — |
| 2 (poly-Q) | QQQQQQQQQQQ | 11 | 7.1e-25 | 4.5 | 0.00022 | **Yes** | — | — |
| 3 (structured) | M[SA]T[TS][IV]METTTT[ML]AT[TS] | 15 | 6.6e-06 | ∞ | 0.086 | No | **98.6%** | DEG_SPOP_SBC_1 |

## Biological Context

- **242 SFARI Category 1 genes** → 1,502 autism IDRs from 217 proteins
- **251 length-matched controls** → 1,130 control IDRs from 196 proteins
- **Three negative controls pass:** shuffled sequences (0 motifs found), positive spike-in (4/4 recovered), bootstrap (Motif 3 stable)
- **Primary finding:** Poly-Q and poly-H tracts are enriched in autism IDRs, likely reflecting the over-representation of transcriptional regulators and chromatin modifiers among high-confidence autism genes
- **Best SLiM candidate:** Motif 3 (SPOP-binding degron-like), which is stable under resampling but underpowered for FDR significance at current sample size
