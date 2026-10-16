from Bio import SeqIO
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

record = SeqIO.read("data/NC_012920.1.gb", "genbank")
seq = str(record.seq).upper()

# GC-skew for Origin OH
window = 500
gc_skew = []
for i in range(0, len(seq)-window, window):
    s = seq[i:i+window]
    g = s.count('G')
    c = s.count('C')
    skew = (g-c)/(g+c) if (g+c)>0 else 0
    gc_skew.append(skew)

min_idx = np.argmin(gc_skew) * window
print(f"OH origin ~ {min_idx} bp via GC-skew minimum")

# Circular map
fig, ax = plt.subplots(figsize=(8,8), subplot_kw={'projection':'polar'})
ax.set_theta_zero_location('N')
for feat in record.features:
    if feat.type == "CDS":
        start = int(feat.location.start)/len(seq)*2*np.pi
        end = int(feat.location.end)/len(seq)*2*np.pi
        ax.barh(1, width=end-start, left=start, height=0.3, color='#2ca02c', alpha=0.8)
    elif feat.type == "tRNA":
        start = int(feat.location.start)/len(seq)*2*np.pi
        end = int(feat.location.end)/len(seq)*2*np.pi
        ax.barh(1, width=end-start, left=start, height=0.3, color='#ff7f0e', alpha=0.8)

ax.set_yticklabels([])
ax.set_title("Human Mito Circular Transcriptome (NC_012920.1)\nGreen=CDS mRNA, Orange=tRNA", pad=20)
Path("figures").mkdir(exist_ok=True)
plt.savefig("figures/circular_map.png", dpi=300, bbox_inches='tight')
print("Saved -> figures/circular_map.png")

# GC-skew plot
plt.figure(figsize=(10,4))
plt.plot(range(0, len(gc_skew)*window, window), gc_skew)
plt.axvline(min_idx, color='red', linestyle='--', label=f'OH ~{min_idx}bp')
plt.xlabel("Position (bp)"); plt.ylabel("GC-skew")
plt.title("GC-skew Origin Detection (OH)")
plt.legend()
plt.savefig("figures/gc_skew_origin.png", dpi=300, bbox_inches='tight')
print("Saved -> figures/gc_skew_origin.png")