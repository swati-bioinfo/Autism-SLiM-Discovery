#!/usr/bin/env python3
"""
Phase 9: Generate publication figures.
- Figure 1: IDR length distribution (autism vs control)
- Figure 2: Motif enrichment bar plot (odds ratios)
- Figure 3: Bootstrap distribution (Motif 3 site count)
"""

import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path

sns.set_theme(style='whitegrid', palette='muted', font_scale=1.1)
plt.rcParams['figure.dpi'] = 150
FIG_DIR = Path('figures')
FIG_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Figure 1: IDR Length Distribution
# ---------------------------------------------------------------------------
def read_idr_lengths(fasta_path):
    lengths = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                lengths.append(len(line.strip()))
    return np.array(lengths)

autism_lengths = read_idr_lengths('00-raw/autism-idrs-clean.fasta')
control_lengths = read_idr_lengths('01-control-data/control-idrs-clean.fasta')

fig, ax = plt.subplots(figsize=(8, 5))
bins = np.logspace(np.log10(15), np.log10(max(autism_lengths.max(), control_lengths.max())), 40)
ax.hist(autism_lengths, bins=bins, alpha=0.6, label=f'Autism (n={len(autism_lengths)})', color='#9F2F2D')
ax.hist(control_lengths, bins=bins, alpha=0.6, label=f'Control (n={len(control_lengths)})', color='#1F6C9F')
ax.set_xscale('log')
ax.set_xlabel('IDR Length (AA, log scale)')
ax.set_ylabel('Count')
ax.set_title('IDR Length Distribution: Autism vs Control')
ax.legend()
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig1-idr-length-distribution.png')
plt.close(fig)
print('Fig 1: IDR length distribution -> figures/fig1-idr-length-distribution.png')

# Stats
print(f'  Autism:  mean={autism_lengths.mean():.1f}  median={np.median(autism_lengths):.1f}  n={len(autism_lengths)}')
print(f'  Control: mean={control_lengths.mean():.1f}  median={np.median(control_lengths):.1f}  n={len(control_lengths)}')

# ---------------------------------------------------------------------------
# Figure 2: Motif Enrichment (Odds Ratios)
# ---------------------------------------------------------------------------
enrich = pd.read_csv('03-analysis/enrichment-results.csv')
# Add a column for plotting
enrich['log_OR'] = np.log2(enrich['odds_ratio'].clip(upper=100))  # cap inf for viz
enrich['sig'] = enrich['q_value'] < 0.05
enrich['neg_log_q'] = -np.log10(enrich['q_value'].clip(lower=1e-10))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Odds ratio bar
colors = ['#1F6C9F' if s else '#C0C0C0' for s in enrich['sig']]
bars = ax1.bar(range(len(enrich)), enrich['log_OR'], color=colors, edgecolor='white', linewidth=0.5)
ax1.set_xticks(range(len(enrich)))
ax1.set_xticklabels([f'M{i+1}' for i in range(len(enrich))])
ax1.set_ylabel('log₂(Odds Ratio)')
ax1.set_title('Motif Enrichment (Fisher\'s Exact Test)')
ax1.axhline(0, color='black', linewidth=0.5)
# Add significance stars
for i, s in enumerate(enrich['sig']):
    if s:
        ax1.text(i, enrich['log_OR'].iloc[i] + 0.1, '*', ha='center', fontsize=16, fontweight='bold', color='#9F2F2D')

# Panel B: Volcano
ax2.scatter(enrich['log_OR'], enrich['neg_log_q'], c=['#9F2F2D' if s else '#C0C0C0' for s in enrich['sig']], s=80)
for i in range(len(enrich)):
    ax2.annotate(f'M{i+1}', (enrich['log_OR'].iloc[i], enrich['neg_log_q'].iloc[i]),
                 textcoords='offset points', xytext=(5, 5), fontsize=8)
ax2.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='q=0.05')
ax2.set_xlabel('log₂(Odds Ratio)')
ax2.set_ylabel('−log₁₀(q-value)')
ax2.set_title('Motif Enrichment Volcano')
ax2.legend()
fig.tight_layout()
fig.savefig(FIG_DIR / 'fig2-motif-enrichment.png')
plt.close(fig)
print('Fig 2: Motif enrichment -> figures/fig2-motif-enrichment.png')
print(enrich[['motif_name', 'odds_ratio', 'p_value', 'q_value', 'sig']].to_string(index=False))

# ---------------------------------------------------------------------------
# Figure 3: Bootstrap Distribution (Motif 3 site count)
# ---------------------------------------------------------------------------
if os.path.exists('03-analysis/bootstrap-results.csv'):
    # Read only the main data rows (before blank line + summary)
    with open('03-analysis/bootstrap-results.csv') as f:
        raw = f.read()
    boot_data = raw.split('\n\n')[0]
    import io
    boot = pd.read_csv(io.StringIO(boot_data))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(boot['sites_found'], bins=range(0, boot['sites_found'].max()+2), alpha=0.7, color='#956400', edgecolor='white')
    ax.axvline(6, color='#9F2F2D', linestyle='--', linewidth=1.5, label='Full set (6 sites)')
    ax.axvline(boot['sites_found'].mean(), color='#1F6C9F', linestyle=':', linewidth=1.5, label=f'Mean ({boot["sites_found"].mean():.1f})')
    ax.set_xlabel('Motif 3 Sites per Bootstrap')
    ax.set_ylabel('Count (out of 1000)')
    ax.set_title('Bootstrap Stability: Motif 3 Site Count Distribution')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'fig3-bootstrap-distribution.png')
    plt.close(fig)
    print(f'Fig 3: Bootstrap distribution -> figures/fig3-bootstrap-distribution.png')
    print(f'  Mean={boot["sites_found"].mean():.1f}, SD={boot["sites_found"].std():.1f}, >=3 sites={boot["sites_found"].ge(3).mean()*100:.1f}%')
else:
    print('Fig 3: bootstrap-results.csv not found, skipping')

print('\nAll figures generated in figures/')
