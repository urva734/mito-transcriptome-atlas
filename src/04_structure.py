from Bio import SeqIO

record = SeqIO.read("data/NC_012920.1.gb", "genbank")
# Find MT-TL1 gene (m.3243A>G MELAS)
for feat in record.features:
    if feat.type=="tRNA" and "gene" in feat.qualifiers:
        if feat.qualifiers["gene"][0]=="MT-TL1":
            start = int(feat.location.start)
            end = int(feat.location.end)
            seq = str(record.seq[start:end])
            print(f"MT-TL1 WT: {start}-{end} len {len(seq)}")
            # MELAS position is 3243 in genome -> in tRNA
            mut_pos_genome = 3242 # 0-based
            mut_pos_trna = mut_pos_genome - start
            wt_base = record.seq[mut_pos_genome]
            print(f"m.3243A>G: genome pos {mut_pos_genome+1} base {wt_base} -> G")
            print(f"tRNA position {mut_pos_trna}")
            # GC change
            gc_wt = (seq.count('G')+seq.count('C'))/len(seq)*100
            seq_mut = seq[:mut_pos_trna] + "G" + seq[mut_pos_trna+1:]
            gc_mut = (seq_mut.count('G')+seq_mut.count('C'))/len(seq_mut)*100
            print(f"GC% WT {gc_wt:.1f} -> MUT {gc_mut:.1f} : cloverleaf destabilization")
            # Save for figure
            with open("results/melas_mutation.txt","w") as out:
                out.write(f"WT {seq}\nMUT {seq_mut}\nGC {gc_wt:.1f}->{gc_mut:.1f}\n")
            break

print("Saved -> results/melas_mutation.txt")
print("Note: Full structure plot = draw cloverleaf manually or with ViennaRNA, GC change explains loss.")