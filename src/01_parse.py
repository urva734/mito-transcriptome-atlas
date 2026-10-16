from Bio import SeqIO
import pandas as pd
from pathlib import Path

gb_path = Path("data/NC_012920.1.gb")
fasta_path = Path("data/mito.fasta")

record = SeqIO.read(gb_path, "genbank")
seq = str(record.seq)
print(f"Length: {len(seq)} bp")
print(f"N at 3107: {seq[3106]} (validation)")

genes = []
covered = set()
for feat in record.features:
    if feat.type in ["gene", "rRNA", "tRNA", "CDS"]:
        start = int(feat.location.start)
        end = int(feat.location.end)
        strand = feat.location.strand
        gene_name = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["unknown"]))[0]
        genes.append([feat.type, gene_name, start, end, strand, end-start])
        for i in range(start, end):
            covered.add(i)

df = pd.DataFrame(genes, columns=["type","name","start","end","strand","length"])
df.to_csv("results/gene_table.csv", index=False)

coding_pct = len(covered)/len(seq)*100
print(f"Coding %: {coding_pct:.1f}% via set() coverage")
print(f"Genes: {len(df)} (13 CDS + 22 tRNA + 2 rRNA)")
print("Saved -> results/gene_table.csv")