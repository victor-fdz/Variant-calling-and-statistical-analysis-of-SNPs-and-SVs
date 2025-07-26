# 🧬 Bioinformatics Internship Scripts

Welcome to this repository!  
Here you'll find some of the scripts I developed during my internship as a **Bioinformatics Research Support** in the **Animal Genomics** area at the **Centre for Research in Agricultural Genomics (CRAG)**.  

---

## 🚀 Projects Overview

### **Project 1: Statistical Detection of GEBV-Correlated Variants**
Scripts for identifying variants correlated with **Genomic Breeding Value (GEBV)**.

- **`Random_Subsets_code.py`**  
  🎲 Random selection of variants to be plotted.

- **`SNPs_Tests_code.py`**  
  📊 Statistical testing of each variant in a VCF, including:  
  - Spearman correlation  
  - Fisher Exact Test  
  - Pearson correlation  
  - Linear Model
    
- **`compare_VCFs_code.py`**  
  🔍 Takes **2 VCF files** and compares their variant content (rows), visualizing the result with a **Venn Diagram**.

---

### **Project 2: Variant Calling & Azoospermia-Causal SV Detection**
Scripts for detecting possible **Structural Variants (SVs)** related to **azoospermia**.

- **`calculate_allele_freq.py`**  
  🧪 Calculates the **allelic frequency** of each variant in a VCF using the number of reads mapped to each option (**var/ref**).

- **`R-code`**  
  🧬 R Scripts and complementary files used to plot a 5-side Venn Diagram of shared variants among animal pools and variants in karyotypes. 
