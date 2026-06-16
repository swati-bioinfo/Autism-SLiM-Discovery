"""Extract SFARI Category 1 (high confidence) gene symbols."""
import csv

input_file = "00-raw/SFARI-Gene_genes.csv"
output_file = "00-raw/category1-genes.txt"

genes = []
with open(input_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("gene-symbol") and row.get("gene-score") == "1":
            genes.append(row["gene-symbol"])

with open(output_file, "w") as f:
    for g in sorted(genes):
        f.write(g + "\n")

print(f"Found {len(genes)} Category 1 genes")
print(f"Saved to {output_file}")
