import base64, os, io, sys
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = 'figures'
os.makedirs(FIG_DIR, exist_ok=True)

tree = ET.parse('03-analysis/meme.xml')
root = tree.getroot()
motifs = root.find('motifs').findall('motif')

for mi, motif in enumerate(motifs):
    e_val = float(motif.get('e_value', '1'))
    name = motif.get('name', f'Motif {mi+1}')
    n_sites = motif.get('sites', '?')
    width = int(motif.get('width', '0'))
    
    print(f'Motif {mi+1}: name="{name}" w={width} sites={n_sites} e={e_val:.0e}')
    
    # Extract PWM
    alphabet_matrix = motif.find('probabilities/alphabet_matrix')
    pwm = []
    for arr in alphabet_matrix.findall('alphabet_array'):
        pos_freqs = {}
        for v in arr.findall('value'):
            letter = v.get('letter_id')
            val = float(v.text)
            if letter and val > 0:
                pos_freqs[letter] = val
        pwm.append(pos_freqs)
    
    # Create sequence logo
    fig, ax = plt.subplots(figsize=(max(width * 0.35, 2.0), 1.5))
    
    # Colors for amino acids (simple scheme)
    aa_colors = {
        'A': '#1f77b4', 'C': '#ff7f0e', 'D': '#e31a1c', 'E': '#d62728',
        'F': '#9467bd', 'G': '#8c564b', 'H': '#e377c2', 'I': '#7f7f7f',
        'K': '#bcbd22', 'L': '#17becf', 'M': '#aec7e8', 'N': '#ffbb78',
        'P': '#98df8a', 'Q': '#ff9896', 'R': '#c5b0d5', 'S': '#c49c94',
        'T': '#f7b6d2', 'V': '#c7c7c7', 'W': '#dbdb8d', 'Y': '#9edae5'
    }
    
    for i in range(width):
        pos = pwm[i]
        # Sort by frequency ascending
        sorted_letters = sorted(pos.keys(), key=lambda l: pos[l])
        bottom = 0
        for letter in sorted_letters:
            freq = pos[letter]
            color = aa_colors.get(letter, '#cccccc')
            ax.bar(i, freq, bottom=bottom, width=0.85,
                   color=color, edgecolor='white', linewidth=0.3)
            if freq > 0.12:
                # Adjust font size based on position density
                fs = max(4, min(7, 7 * freq))
                ax.text(i, bottom + freq/2, letter, ha='center', va='center',
                        fontsize=fs, fontweight='bold', color='white')
            bottom += freq
    
    ax.set_xticks(range(width))
    ax.set_xticklabels(range(1, width+1), fontsize=5)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Prob', fontsize=6)
    e_str = f'{e_val:.1e}'
    ax.set_title(f'Motif {mi+1} (E={e_str})', fontsize=7)
    if width > 20:
        ax.xaxis.set_tick_params(labelsize=3)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    fig.tight_layout(pad=0.3)
    figpath = os.path.join(FIG_DIR, f'motif{mi+1}-logo.png')
    fig.savefig(figpath, dpi=150)
    plt.close(fig)
    print(f'  -> saved {figpath}')
