"""Create positive control by spiking known SLiMs into synthetic IDR sequences.

Adds sequences containing known short linear motifs to the autism IDR set
to validate that the MEME pipeline can detect biologically meaningful motifs.

Known SLiMs spiked:
  - SH3-binding Class I:   PxxPxR/K  (ELM: LIG_SH3_1)
  - SH3-binding Class II:  R/KxxPxxP (ELM: LIG_SH3_2)
  - 14-3-3 binding:        R[S/Ar][S/Ar][S/Ar]P (ELM: LIG_14-3-3_1)
  - PDZ-binding (Class I): x[S/T]x[L/V]-COOH (ELM: LIG_PDZ_1)

Each spiked sequence: random IDR-like background (P,S,Q,E-rich) + embedded motif.
Total: 100 positive-control sequences added to autism set, 100 to control set.

Output: 00-raw/autism-idrs-positive-control.fasta
        01-control-data/control-idrs-positive-control.fasta
"""
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

random.seed(7)

def idr_like_background(length=60):
    """Generate IDR-like background sequence (P,S,Q,E,T,A rich)."""
    # Frequencies approximating IDR composition
    freqs = {
        'P': 0.12, 'S': 0.12, 'Q': 0.08, 'E': 0.08, 'T': 0.07,
        'A': 0.07, 'G': 0.07, 'K': 0.06, 'R': 0.05, 'D': 0.05,
        'N': 0.04, 'L': 0.04, 'V': 0.03, 'H': 0.03, 'Y': 0.02,
        'F': 0.02, 'M': 0.02, 'I': 0.02, 'C': 0.01, 'W': 0.01,
    }
    aas = list(freqs.keys())
    weights = list(freqs.values())
    return ''.join(random.choices(aas, weights=weights, k=length))

motifs = {
    'SH3_ClassI': {
        'pattern': 'PxxPxR',
        'consensus': 'PxxPxR',
        'notes': 'SH3-binding Class I (ELM:LIG_SH3_1)',
    },
    'SH3_ClassII': {
        'pattern': 'RxxPxxP',
        'consensus': 'RxxPxxP',
        'notes': 'SH3-binding Class II (ELM:LIG_SH3_2)',
    },
    '14-3-3': {
        'pattern': 'RSxSP',
        'consensus': 'RSxSP',
        'notes': '14-3-3 binding mode 1 (ELM:LIG_14-3-3_1)',
    },
    'PDZ_ClassI': {
        'pattern': 'STxL',
        'consensus': 'SxTL',
        'notes': 'PDZ-binding Class I C-term (ELM:LIG_PDZ_1)',
    },
}

# For each motif: generate spiked sequences
# Half embedded in the middle, half at varying positions
spiked_autism = []
spiked_control = []

for motif_name, motif_info in motifs.items():
    pattern = motif_info['pattern']
    for i in range(25):  # 25 per motif = 100 total
        bg = idr_like_background(random.randint(40, 80))
        insert_pos = random.randint(5, len(bg) - 10)
        spiked_seq = bg[:insert_pos] + pattern + bg[insert_pos:]
        
        rec = SeqRecord(
            Seq(spiked_seq),
            id=f"POSCTRL_{motif_name}_{i+1}_{pattern}",
            description=f"Positive control: {motif_info['notes']}",
        )
        spiked_autism.append(rec)
        spiked_control.append(rec)

# Write positive control files
from Bio import SeqIO

# For autism: merge spiked with original
orig_autism = list(SeqIO.parse("00-raw/autism-idrs-clean.fasta", "fasta"))
SeqIO.write(orig_autism + spiked_autism,
            "02-motif-discovery/autism-idrs-positive-control.fasta", "fasta")
print(f"Autism: {len(orig_autism)} original + {len(spiked_autism)} spiked = {len(orig_autism)+len(spiked_autism)}")

# For control: merge spiked with original
orig_control = list(SeqIO.parse("01-control-data/control-idrs-clean.fasta", "fasta"))
SeqIO.write(orig_control + spiked_control,
            "02-motif-discovery/control-idrs-positive-control.fasta", "fasta")
print(f"Control: {len(orig_control)} original + {len(spiked_control)} spiked = {len(orig_control)+len(spiked_control)}")

print("\nPositive control sequences created.")
print("To validate: submit to MEME discriminative mode.")
print("Expected: MEME should detect PxxPxR and RxxPxxP motifs.")
