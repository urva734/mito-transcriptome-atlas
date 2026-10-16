# Module 04: tRNA Punctuation Model

Model: Single polycistronic transcript -> endonucleolytic cleavage at tRNAs -> 37 mature RNAs

Evidence: src/03_rna_logic.py output simulates 1 long RNA (~16kb) cleaved at 22 tRNA excision points -> 37 discrete RNAs

Mechanism:
- tRNAs act as punctuation marks via cloverleaf structure
- RNase P cleaves 5' end, RNase Z (ELAC2) cleaves 3' end
- Releases flanking mRNAs/rRNAs

Conclusion: Explains how overlapping intronless transcripts generate individual RNAs without splicing.