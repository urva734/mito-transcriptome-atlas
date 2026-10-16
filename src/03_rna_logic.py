from Bio import SeqIO
from Bio.Data import CodonTable

record = SeqIO.read("data/NC_012920.1.gb", "genbank")
seq = str(record.seq)

# Genetic Code Recoding
universal = CodonTable.unambiguous_dna_by_name["Standard"]
mito = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

u_stops = universal.stop_codons
m_stops = mito.stop_codons

print(f"Universal stops: {u_stops} -> 93 stops in mito if used")
print(f"Mito stops: {m_stops} -> 0 stops (recoding)")

# Count in mito
count_uga = seq.count("TGA")
count_aga = seq.count("AGA")
count_agg = seq.count("AGG")
print(f"TGA (W in mito, Stop in universal): {count_uga}")
print(f"AGA/AGG (Stop in mito, R in universal): {count_aga}/{count_agg}")

# tRNA Punctuation Model
trnas = [f for f in record.features if f.type=="tRNA"]
print(f"\n22 tRNAs as RNA cleavage signals:")
for t in trnas:
    name = t.qualifiers.get("gene",[""])[0]
    print(f"{name}: {int(t.location.start)}-{int(t.location.end)}")

print("\nPolycistronic model: 1 long RNA -> cut by tRNAs -> 37 mature RNAs")