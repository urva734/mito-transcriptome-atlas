# 🧬 Human Mitochondrial Transcriptome Atlas

An end-to-end pure-Python pipeline that annotates human mitochondrial genome (16,569 bp) from genome architecture to disease mapping. Built from first principles without Biopython.

## 👤 Author
**Urva Sohail**

## ✨ Features - 5 Modules in 1 Repo
- **🧬 Module 01**: Genome architecture + circular map
- **📈 Module 02**: GC-skew origin of replication detection
- **🔄 Module 03**: Genetic recoding UGA->Trp mapping
- **✂️ Module 04**: tRNA punctuation model simulation
- **🧬 Module 05**: MELAS disease m.3243A>G analysis

## 💻 Technologies Used
- **Python 3+**: Core language
- **math + pathlib**: File handling + polar math
- **matplotlib**: For circular_map + gc_skew plots
- **pandas + numpy**: For tables + sliding windows
- **Custom Parsers**: FASTA parser without Biopython

## 🚀 How to Run
Using VS Code:
1. Clone repo `git clone https://github.com/urva734/mito-transcriptome-atlas.git`
2. Install deps `pip install -r requirements.txt`
3. Run workflow:
   ```
python src/01_parse.py
python src/02_visualize.py
python src/03_rna_logic.py
python src/04_structure.py
```
4. Or run module-wise from `workflow/module_01_architecture/` to `module_05_disease/`

## 📥 Sample Input & Output
**Input:** Human mito FASTA 16,569 bp (NC_012920.1)

**Output Module 01:**
```
Length: 16569 bp | Genes: 37 (13 CDS, 22 tRNA, 2 rRNA)
Coding: 92.7% | N at 3107 reference ambiguity
Circular map generated
```
**Output Module 02:**
```
OH origin ~14000 bp D-loop | OL origin ~5800 bp WANCY
Heavy strand G-rich confirmed
```

**Output Module 03-05:**
```
UGA->Trp: 12 | AUA->Met: 45 | m.3243A>G G-U wobble destabilized
1 polycistronic RNA -> 22 tRNA cuts -> 37 mature RNAs
```

## 📊 Generated Plots
`figures/circular_map.png`:
- X = Circular genome 0-16569 bp
- Colors = CDS, tRNA, rRNA - Shows compact 92.7% coding

`figures/gc_skew_origin.png`:
- X = Position, Y = GC-skew (G-C)/(G+C)
- Peak at 14000 bp = OH origin D-loop

## 📁 Downloaded Files
```
figures/circular_map.png - Circular genome map
figures/gc_skew_origin.png - GC-skew origin plot
src/01_parse.py - FASTA parser
src/02_visualize.py - Circular + skew plots
src/03_rna_logic.py - Recoding logic
src/04_structure.py - Punctuation + disease
workflow/module_01 to 05 - Detailed modules
LICENSE - MIT license
```

## ⚙️ How It Works
- **01:** Parses FASTA, extracts 37 genes, polar coordinate circular plot
- **02:** Counts G/C per 500bp window, (G-C)/(G+C), peak detection
- **03:** Compares universal vs mito code, counts recoded codons
- **04:** Simulates 16kb polycistronic cut at 22 tRNAs -> 37 RNAs
- **05:** Maps m.3243A>G, A-U vs G-U stability, OXPHOS impact

## 📊 Result Summary
- Total Length: 16569 bp parsed
- Genes: 37 mapped | Coding Density: 92.7%
- Origins: OH ~14000 bp, OL ~5800 bp
- Recodings: 4 types validated
- Disease: MELAS mechanism explained
- Status: PASS - Complete atlas

## 🔮 Future Improvements
- [ ] Add interactive Plotly circular map
- [ ] Add more disease variants (LHON, MERRF)
- [ ] Add AT-skew analysis
- [ ] Export GenBank format

