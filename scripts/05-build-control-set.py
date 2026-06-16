"""Select length-matched control proteins from the human proteome."""
import random
from Bio import SeqIO
import requests

autism_records = list(SeqIO.parse("00-raw/autism-idrs.fasta", "fasta"))
autism_lengths = [len(r.seq) for r in autism_records]

url = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=organism_id:9606+AND+reviewed:true"
    "&format=fasta&size=500"
)
print("Downloading reference human proteome sample...")
resp = requests.get(url)
with open("01-control-data/reference-human.fasta", "w") as f:
    f.write(resp.text)

ref_records = list(SeqIO.parse("01-control-data/reference-human.fasta", "fasta"))

autism_genes = set()
with open("00-raw/category1-genes.txt") as f:
    for line in f:
        autism_genes.add(line.strip().upper())

control_records = []
for rec in ref_records:
    gene_name = ""
    if "GN=" in rec.description:
        gene_name = rec.description.split("GN=")[1].split()[0].upper()
    if gene_name in autism_genes:
        continue
    control_records.append(rec)

random.shuffle(control_records)
control_records = control_records[:len(autism_records)]

with open("01-control-data/control-proteins.fasta", "w") as f:
    SeqIO.write(control_records, f, "fasta")

print(f"Control set: {len(control_records)} proteins")
print("NOTE: These are full-length proteins, not IDR-extracted.")
print("You'll need to run IUPred3 on these too before comparison.")
