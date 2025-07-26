## 1. Plotting Structural Variants in a Karyotype (📁 `SVs_Plots`)
In this analysis I'm searching 🔎 for candidate SVs that can be causing an specific disease in an specific animal specie. The file `Plot_SVs_Script.R` ⌨️ is an script that takes as arguments: 
- 📄 `SVs.txt`: an example file with randomly generated SVs.
- 📄 `chr_tab_len.txt`: each chromosome of the reference genome with its length. Also artificially modified to keep real data private.
  
The script generates 3 plots, whose examples are also in the folder: 
- 🎨 `plot1.txt`: SVs in the reference genome.
- 🎨 `plot2.txt`: high density 📈 zone of SVs.
- 🎨 `plot3.txt`: high density zone zoomed, with SVs types labeled. 

## 2. Plotting Venn Diagrams of Private/Shared SVs among Animal Genetic Pools (📁 `Venn_Diagrams`)

In this second case I'm comparing animal genetic pools 🧬. After calling the SVs, I made comparisons of 2 pools, but I wanted to do bigger comparisons (between >2 pools at the same time). That's why I created `Venn_Plot.R` ⌨️, an R script that takes as arguments pairs of files: 
- 📄 `SVs.txt`: example random file with SVs that has the number of alleles for 2 pools (initial comparison). 
- 📄 `stats.txt`: example random file with the number of statistically significant SVs in this initial comparison.

The script permits a **naive visualisation of how many SVs are private and how many are shared between genetic pools**. It is actually coded to give 10 arguments (5 pairs), but can be easily modified 🔧 to reduce the number of pools to be compared. 
