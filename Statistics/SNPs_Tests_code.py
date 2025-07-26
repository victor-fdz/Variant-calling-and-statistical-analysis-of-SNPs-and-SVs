# import necessary modules
import numpy as np
import pandas as pd
import scipy.stats as sci
import statsmodels.api as sm
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection
import sys

#-#-#-#-#-#-#-#-#-#-#-# FILE ACCESSION, READS to GENOMES CONVERSION & FISHER/SPEARMAN TESTING #-#-#-#-#-#-#-#-#-#-#-#
# Define files
file_path = sys.argv[1]

# Open it and save its lines
with open(file_path) as file: 
    lines = file.readlines()

# Search the first variant line
for i in range(0,len(lines)): 
    if not lines[i].startswith("#"): # If the line do not start with #, it does not belong to the header
        col_names = lines[i-1]
        z = i
        break

# Define the dataframe
df = pd.DataFrame([i[:-1].split("\t") for i in lines[z:]], # z is the index where the header finishes
                  columns = col_names[1:-1].split("\t")) # First character is a "#" and last one is a line break
        
    
# Define important lists
samples = ["ABA09619","ABA09622","ABA09623","ABA09951","ABA09952","ABA09953","ABA09954","ABA09955"] # All samples' names
samples_12 = ["ABA09623","ABA09951","ABA09952","ABA09953"] # Samples with 12 animals  (24 genomes)
GEBVs = [-4.401016364, -3.184106364, -2.452780833, -2.081470833, 3.170856667, 3.466019167, 3.997353636, 5.002573636] # Genomic breeding values
genomes = [] # Empty list to append the number of genomes for each SNP and each sample
g_freqs = [] # Empty list to store the frequencies of genomes
spearman_results = [] # Empty list to store the results of Spearman Tests
pearson_results = [] # Empty list to store the results of Spearman Tests
reg_results = [] # Empty list to store the results of Spearman Tests
fisher_results = [] # Empty list to store the results of Fisher Exact Tests between superpools
fisher_18_results = [] # Empty list to store the results of Fisher Exact Tests between L1 and H8
fisher_45_results = [] # Empty list to store the results of Fisher Exact Tests between L4 and H5

# Access to each SNP and test it
for i in range(0, len(df)): 
    SNP = df.iloc[i] # Define the SNP
    genomes_SNP = [] # Empty list to store the number of genomes for each SNP
    
    # Empty lists to store the number of genomes of superpools H and L
    ref_L = []
    ref_H = []
    alt_L = []
    alt_H = []
    
    # Access to each sample
    for s in range(0,len(samples)): 
        reads = SNP[samples[s]].split(":")[1].split(",") # Define its reads
        ref = int(reads[0])
        alt = int(reads[1])
        
        # Number of animals per pool
        n = 11 # Default
        if samples[s] in samples_12: # If the sample is in the list of samples with 12 animals
            n = 12 # Change n
        
        # Filter to avoid "div by 0" error
        if (ref + alt) == 0: # If the sample do not have reads for that SNP, define genomes as 0 directly
            g_alt = 0
            g_ref = 0
        
        else: # Convert reads into genomes
            g_ref = round((ref / (ref + alt)) * (n * 2), 0) # n*2 because pigs are diploids
            g_alt = round((alt / (ref + alt)) * (n * 2), 0)
        
        # Append the genomes of the sample to the SNP list
        genomes_SNP.append((g_ref, g_alt)) 
        
        # Add the number of genomes to the corresponding superpool
        if s < 4: # For the first 4 samples (Low superpool)
            ref_L.append(g_ref)
            alt_L.append(g_alt)
        else: 
            ref_H.append(g_ref)
            alt_H.append(g_alt)

    # Add the number of ref/alt genomes of the superpools to the genomes list
    low = (sum(ref_L), sum(alt_L))
    high = (sum(ref_H), sum(alt_H))
    genomes_SNP.extend([low, high])  
        
    # Append the genomes of the SNP to the general genomes list   
    genomes.append(genomes_SNP) 
    
    # Fisher Exact Test Low vs High (superpools)
    g_fisher = sci.fisher_exact(np.array([low, high]))
    fisher_results.append(g_fisher)
    
    # Fisher Exact Test L1 vs H8
    g_fisher = sci.fisher_exact(np.array([genomes_SNP[0], genomes_SNP[-3]]))
    fisher_18_results.append(g_fisher)

    # Fisher Exact Test L4 vs H5
    g_fisher = sci.fisher_exact(np.array([genomes_SNP[3], genomes_SNP[4]]))
    fisher_45_results.append(g_fisher)

    # Spearman Test
    # Transform the numbers of genomes into frequencies 
    freqs = []
    for pair in genomes_SNP[:8]:
        if pair[0] + pair[1] == 0:
            f = 0
        else: 
            f = pair[1] / (pair[0] + pair[1]) # alt / (ref + alt)
        freqs.append(f)
    g_freqs.append(freqs)
    
    # Filter: if the frequency is constant across pools, avoid testing
    if len(set(freqs)) == 1: 
        spearman_results.append((0, 1)) # Correlation = 0; p-value = 1
        pearson_results.append((0, 1)) # Correlation = 0; p-value = 1
        reg_results.append(".")
        continue
    
    # Perform Spearman test
    g_spearman = sci.spearmanr(np.array(GEBVs), np.array(freqs))
    spearman_results.append(g_spearman)
    
    # Pearson test (all pools)
    g_pearson = pearsonr(np.array(GEBVs), np.array(freqs))
    pearson_results.append(g_pearson)
    
    # Linear Model test (all pools)
    if g_pearson[1] < 0.05:
        x_with_constant = sm.add_constant(GEBVs)
        LM = sm.OLS(freqs,x_with_constant).fit()
        reg_results.append([LM.params,LM.pvalues])
    else: 
        reg_results.append(".")
    
        
        
# Changing format for better visualization
new_reg = []
for i in reg_results:
    if i != ".":
        new_reg.append(";".join([",".join([str(a) for a in i[0]]),",".join([str(a) for a in i[1]])]))
    else:
        new_reg.append(".")
    
genomes_col = [] # Empty list to add the genomes better formatted for visualization

# Format transformation
for i in range(0,len(genomes)): 
    string = ""
    for s in range(0,len(genomes[i])):
        string += str(str(int(genomes[i][s][0]))+","+str(int(genomes[i][s][1])))
        if not s == len(genomes[i])-1:
            string += ":"
    genomes_col.append(string)
    
# Add the new columns to the dataframe
df["Genomes_Samples&Superpools"] = genomes_col
df["Frequencies(Genomes)"] = [":".join([str(a) for a in i]) for i in g_freqs]
    
# FDR correction
reject_f, corrected_f = fdrcorrection([float(i[1]) for i in fisher_results], alpha = 0.05)
reject_f18, corrected_f18 = fdrcorrection([float(i[1]) for i in fisher_18_results], alpha = 0.05)
reject_f45, corrected_f45 = fdrcorrection([float(i[1]) for i in fisher_45_results], alpha = 0.05)
reject_s, corrected_s = fdrcorrection([float(i[1]) for i in spearman_results], alpha = 0.05)
reject_p, corrected_p = fdrcorrection([float(i[1]) for i in pearson_results], alpha = 0.05)

# Add tests results
df["Spearman_Correlation"] = [i[0] for i in spearman_results]
df["Spearman_qval"] = corrected_s

df["Fisher_LvH_OddsRatio"] = [i[0] for i in fisher_results]
df["Fisher_LvH_qval"] = corrected_f

df["Fisher_L1vH8_OddsRatio"] = [i[0] for i in fisher_18_results]
df["Fisher_L1vH8_qval"] = corrected_f18

df["Fisher_L4vH5_OddsRatio"] = [i[0] for i in fisher_45_results]
df["Fisher_L4vH5_qval"] = corrected_f45

df["Pearson_r"] = [i[0] for i in pearson_results]
df["Pearson_qvalue"] = corrected_p

df["LinearModel"] = new_reg

# Show the final table
df

# Save it into a file
with open("SNPs_Tests_results.txt", 'w') as f:
    cols = "\t".join(list(df.columns))
    f.write(cols+"\n")
    for i in range(0,len(df)):
        f.write("\t".join([str(a) for a in df.iloc[i]])+"\n")
