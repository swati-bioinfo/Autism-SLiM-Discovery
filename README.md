# Autism Motif Discovery Project

## Overview

Systematic discovery of conserved short linear motifs (SLiMs) in intrinsically disordered
regions (IDRs) of high-confidence autism risk genes (SFARI Category 1, n=242) versus
length-matched controls (n=251 non-SFARI Swiss-Prot human proteins).

## Pipeline

**1**. **Data Acquisition** - SFARI Gene List + UniProt Swiss-Prot sequences
**2**. **Quality Control** - Sequence validation, length distribution, AA composition
**3**. **IDR Prediction** - IUPred3 REST API (score > 0.5, min 15 AA)
**4**. **Control Set** -  Length-matched non-SFARI proteins (±10–50%, seed=42)
**5**. **Motif Discovery** - MEME v5.5.9 discriminative mode (case/control)
**6**. **Validation** - TOMTOM vs ELM, bootstrap (1000×), shuffled control, spike-in control
**7**. **Enrichment** - Fisher's exact test + Benjamini-Hochberg FDR

## Key Results

| Result | Detail |
|---|---|
| Autism IDRs | 1,502 from 217/242 proteins |
| Control IDRs | 1,130 from 196/251 proteins |
| Significant motifs | 3 (poly-H, poly-Q, structured M-S/T-T/S-I/V-M-E-T-T-T-T-M/L-A-T/S) |
| Validated stable | Motif 3 bootstrap 98.6% detection (1000×) |
| TOMTOM best hit | DEG_SPOP_SBC_1 (SPOP-binding degron, p=0.005) |
| Shuffled control | 0 significant motifs (not due to composition noise) |
| Positive control | 4/4 spiked SLiMs recovered (pipeline validates) |
| **Enriched (FDR q<0.05)** | poly-Q: OR=4.5, q=0.00022; poly-H: OR=∞, q=0.0059 |

## Reproducibility

### Data Sources
- SFARI Gene List: https://gene.sfari.org/human-gene (`00-raw/SFARI-Gene_genes.csv`)
- UniProt Swiss-Prot: `rest.uniprot.org/uniprotkb/stream` (`00-raw/uniprot_sprot_human.fasta`)
- IUPred3: `iupred3.elte.hu/iupred3/{accession}`

### Dependencies
- Python 3.10+ with: requests, biopython, pandas, numpy, scipy
- MEME Suite v5.5.9: https://meme-suite.org/meme/tools/meme
- TOMTOM: https://meme-suite.org/meme/tools/tomtom

### Random Seeds
- Control selection: 42
- Shuffle control: 99
- Spike-in control: 7
- Bootstrap: 42

## Directory Structure

```
00-raw/                  # Raw data: SFARI CSV, UniProt sequences, IUPred3 IDRs
01-control-data/         # Length-matched controls, shuffled controls
02-motif-discovery/      # MEME inputs (positive control spike-in)
03-analysis/             # MEME XML, TOMTOM results, bootstrap, enrichment
04-docs/                 # Documentation: index.html, project-steps.md, parameters.log, qc-summary.txt
scripts/                 # All analysis scripts (00-*.py through 09-*.py)
```

## Full Report

Open **[interactive pipeline report](https://swatibio.github.io/Autism-SLiM-Discovery/)** in a browser.

## License

This work is licensed under CC BY-NC-ND 4.0. See `LICENSE` for details.
Copyright (c) 2026 Swati. All rights reserved.
