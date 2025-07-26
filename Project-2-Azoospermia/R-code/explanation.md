# "Calaix de Sastre" - Víctor Fernández 💻
In this repository you will find cases where I used what we have learnt in the practical cases (**ggplot2** R tool ⚙️, **numpy** and **pandas** Python packages...) in a real data analysis that I'm currently doing in **Center for Research on agricultural Genomics** (CRAG), where I'm doing my second internship 🔬.

## 1. Plotting Structural Variants in a Karyotype (📁 `SVs_Plots`)
In this analysis I'm searching 🔎 for candidate SVs that can be causing an specific disease in an specific animal specie. The file `Plot_SVs_Script.R` ⌨️ is an script that takes as arguments: 
- 📄 `SVs.txt`: an example file with randomly generated SVs.
- 📄 `chr_tab_len.txt`: each chromosome of the reference genome with its length. Also artificially modified to keep real data private.

To obtain these files, I used bash coding. First of all, I downloaded the reference genome of my interested specie from NCBI (*genome_info.tsv*). I already had the SVs called (*Candidate_SVs.txt*).
- **head -n 21 genome_info.tsv | tail -n 20 | cut -f3,11 > chr_tab_len.txt** / To take only the chromosomes of the reference genome (not the unmapped contigs) and only the name/number of the chromosome (3rd column) and its length (11th column).
- **head -n 52 Candidate_SVs.txt | tail -n 50 |  cut -f1,2,3,9 > SVs.txt** / To ignore the heading and only taking the chromosome where the SV is (1st column), its start coord (2nd), its length (3rd) and its type (9th).

The script generates 3 plots, whose examples are also in the folder: 
- 🎨 `plot1.txt`: SVs in the reference genome.
- 🎨 `plot2.txt`: high density 📈 zone of SVs.
- 🎨 `plot3.txt`: high density zone zoomed, with SVs types labeled. 

## 2. Plotting Venn Diagrams of Private/Shared SVs among Animal Genetic Pools (📁 `Venn_Diagrams`)

In this second case I'm comparing animal genetic pools 🧬. After calling the SVs, I made comparisons of 2 pools, but I wanted to do bigger comparisons (between >2 pools at the same time). That's why I created `Venn_Plot.R` ⌨️, an R script that takes as arguments pairs of files: 
- 📄 `SVs.txt`: example random file with SVs that has the number of alleles for 2 pools (initial comparison). 
- 📄 `stats.txt`: example random file with the number of statistically significant SVs in this initial comparison.

The script permits a **naive visualisation of how many SVs are private and how many are shared between genetic pools**. It is actually coded to give 10 arguments (5 pairs), but can be easily modified 🔧 to reduce the number of pools to be compared. 

A plot example is also available (🎨 `Venn_Diagram.txt`). 

Although I didn't use **ggplot2** here, the plotting package had a similar usage, so having previous knowledge was really useful.

## 3. Using Numpy/Pandas to Store/Visualize VCF data  (📁 `np_pd_VCF`)

In this final case I wanted to better visualize my VCF files, so I applied these 2 packages that we used in the last practical sessions. The example of the session was really useful 💡, because I had never used this type of packages before in Python. 

The script `VCF_df.py` ⌨️ takes a VCF file, which in this case is the same randomly generated as in the first analysis (📄 `SVs.txt`) and converts its data into a dataframe. How is visualized in a Jupyter Notebook is shown in 🎨 `screenshot_df.png`.
