"""Fetch reviewed (Swiss-Prot) sequences for a list of human genes."""
import requests
import time

with open("00-raw/category1-genes.txt") as f:
    genes = [line.strip() for line in f if line.strip()]

output_file = "00-raw/autism-proteins.fasta"
skipped = []

with open(output_file, "w") as out:
    for i, gene in enumerate(genes):
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query=gene:{gene}+AND+organism_id:9606+AND+reviewed:true"
            "&format=fasta"
        )
        resp = requests.get(url)
        if resp.status_code == 200 and resp.text.strip():
            out.write(resp.text)
            if not resp.text.startswith(">"):
                skipped.append(gene)
        else:
            skipped.append(gene)
        time.sleep(0.2)
        if (i + 1) % 50 == 0:
            print(f"Processed {i+1}/{len(genes)} genes")

if skipped:
    print(f"Warning: {len(skipped)} genes had no Swiss-Prot entry:")
    for g in skipped[:10]:
        print(f"  {g}")
    if len(skipped) > 10:
        print(f"  ... and {len(skipped)-10} more")
print(f"Sequences saved to {output_file}")
