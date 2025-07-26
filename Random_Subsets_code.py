# Import necessary modules
import numpy as np
import pandas as pd
import scipy.stats as sci
from statsmodels.stats.multitest import fdrcorrection
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from matplotlib.pyplot import savefig
import random
plt.style.use('bmh')
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
df = pd.DataFrame([i.split("\t") for i in lines[z:]], # z is the index where the header finishes
                  columns = col_names[1:-1].split("\t")) # First character is a "#" and last one is a line break
        
    
# Define important lists
samples = ["ABA09619","ABA09622","ABA09623","ABA09951","ABA09952","ABA09953","ABA09954","ABA09955"] # All samples' names
samples_12 = ["ABA09623","ABA09951","ABA09952","ABA09953"] # Samples with 12 animals  (24 genomes)
GEBVs = [-4.401016364, -3.184106364, -2.452780833, -2.081470833, 3.170856667, 3.466019167, 3.997353636, 5.002573636] # Genomic breeding values
genomes = [] # Empty list to append the number of genomes for each SNP and each sample
g_freqs = [] # Empty list to store the frequencies of genomes
fisher_results = [] # Empty list to store the results of Fisher Exact Tests between superpools
fisher_18_results = [] # Empty list to store the results of Fisher Exact Tests between L1 and H8
fisher_45_results = [] # Empty list to store the results of Fisher Exact Tests between L4 and H5
spearman_results = [] # Empty list to store the results of Spearman Tests

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
        continue
    
    # Perform Spearman test
    g_spearman = sci.spearmanr(np.array(GEBVs), np.array(freqs))
    spearman_results.append(g_spearman)
    
genomes_col = [] # Empty list to add the genomes better formatted for visualization

# Format transformation
for i in range(0,len(genomes)): 
    string = ""
    for s in range(0,len(genomes[i])):
        string += str(str(int(genomes[i][s][0])) + ","+str(int(genomes[i][s][1])))
        if not s == len(genomes[i])-1:
            string += ":"
    genomes_col.append(string)
    
# Add the new column to the dataframe
df["Genomes_Samples&Superpools"] = genomes_col
df["Frequencies(Genomes)"] = [":".join([str(a) for a in i]) for i in g_freqs]
    
# FDR correction
reject_f, corrected_f = fdrcorrection([float(i[1]) for i in fisher_results], alpha = 0.05)
reject_f18, corrected_f18 = fdrcorrection([float(i[1]) for i in fisher_18_results], alpha = 0.05)
reject_f45, corrected_f45 = fdrcorrection([float(i[1]) for i in fisher_45_results], alpha = 0.05)
reject_s, corrected_s = fdrcorrection([float(i[1]) for i in spearman_results], alpha = 0.05)


# Add Spearman & Fisher tests' results
df["Spearman_Correlation"] = [i[0] for i in spearman_results]
df["Spearman_pval"] = [i[1] for i in spearman_results]
df["Spearman_qval"] = corrected_s
df["Spearman_rejectNull"] = reject_s
df["Fisher_OddsRatio"] = [i[0] for i in fisher_results]
df["Fisher_pval"] = [i[1] for i in fisher_results]
df["Fisher_qval"] = corrected_f
df["Fisher_rejectNull"] = reject_f
df["Fisher_L1vH8_qval"] = corrected_f18
df["Fisher_L4vH5_qval"] = corrected_f45
    
#-#-#-#-#-#-#-#-#-#-#-# FILTERING for SIGNIFICANT VARIANTS and VENN DIAGRAM plot #-#-#-#-#-#-#-#-#-#-#-#
    
# Filter for significance in Spearman Test, Fisher Exact Test and both 
sig_spearman = df[df["Spearman_rejectNull"] == True]
sig_fisher = df[df["Fisher_rejectNull"] == True]
sig_both = sig_fisher[sig_fisher["Spearman_rejectNull"] == True]

# Generate an Excel with significant variants
sig_both.to_excel("SNPs_SignificantSpearman&Fisher.xlsx", index = False)

# Plot Venn Diagram
fig, ax = plt.subplots()
sp = len(sig_spearman)-len(sig_both) # N variants significative only for Spearman
fish = len(sig_fisher)-len(sig_both) # N variants significative only for Fisher
both = len(sig_both) # N variants significative for both tests
ven = venn2(subsets = (sp, fish, both), set_labels = ('Spearman', 'Fisher'))

# Do not plot the numbers
for idx in ('10', '01', '11'):
    if ven.get_label_by_id(idx):
        ven.get_label_by_id(idx).set_text('')

plt.title("Significant SNPs")
plt.text(0.1 ,0.1, f"Only Spearman: {sp}.\nOnly Fisher: {fish}.\nBoth: {both}.")
savefig("Significant_SNPs.png")

# Second plot: we need new subsets
sig_f18 = df[df["Fisher_L1vH8_qval"] < 0.05]
sig_f45 = df[df["Fisher_L4vH5_qval"] < 0.05]
sig_f_f18 = sig_fisher[sig_fisher["Fisher_L1vH8_qval"] < 0.05]
sig_f_f45 = sig_fisher[sig_fisher["Fisher_L4vH5_qval"] < 0.05]
sig_f18_f45 = sig_f18[sig_f18["Fisher_L4vH5_qval"] < 0.05]
sig_all = sig_f18_f45[sig_f18_f45["Fisher_rejectNull"] == False]

# Plot Venn Diagram (Fishers)
fig, ax = plt.subplots()
# A = L1vH8
# B = L4vH5
# C = LvH
ABc = len(sig_f18_f45) - len(sig_all)
AbC = len(sig_f_f18) - len(sig_all)
aBC = len(sig_f_f45) - len(sig_all)
ABC = len(sig_all)
Abc = len(sig_f18) - ABc - AbC - ABC
aBc = len(sig_f45) - ABc - aBC - ABC
abC = len(sig_fisher) - aBC - AbC - ABC

ven = venn3(subsets = (Abc, aBc, ABc, abC, AbC, aBC, ABC), set_labels = ('FisherL1vH8', 'FisherL4vH5','FisherLvH'))
for idx in ('100', '110', '111','101', '001', '011','010'):
    if ven.get_label_by_id(idx):
        ven.get_label_by_id(idx).set_text('')

plt.title("Significant SNPs")
plt.text(0.6,0.5,f"Only LvH: {abC}.\nOnly L1vH8: {Abc}.\nOnly L4vH5: {aBc}.\nAll: {ABC}.")
savefig("venn_fishers.png")

#-#-#-#-#-#-#-#-#-#-#-# DISTRIBUTION and CORRELATION between Q-VALUES #-#-#-#-#-#-#-#-#-#-#-#

# Define the figure
fig, ax = plt.subplots(2,2,figsize=(12,8))

# Extract q-values into lists
s_q = list(df["Spearman_qval"]) # Spearman
f_q = list(df["Fisher_qval"]) # Fisher 

# Define the statistical threshold
threshold  = 0.05

# q-VALUE FISHER DISTRIBUTION          

# Define an histogram and extract its values
n, bins, rect = ax[0,0].hist(f_q, bins = 30)

# Assign colors based on bin values
for rect, bin_val in zip(rect, bins): # For each rectangle (patch) and value of this rectangle
    if bin_val < threshold: # If the value is lower than the threshold
        rect.set_facecolor('forestgreen') # Significant
        rect.set_label("Significant")
    else: # Otherwise
        rect.set_facecolor('grey')  # No
        rect.set_label("Non-significant")

# Axes settings
ax[0,0].set_title("Fisher q-values Distribution", fontsize = 10)
ax[0,0].set_xlabel("Fisher q-val", fontsize = 10)
ax[0,0].set_ylabel("Frequency", fontsize = 10)
ax[0,0].tick_params(labelsize = 8)
ax[0,0].axvline(0.05, color = 'gray', linewidth = 1, linestyle = '--', label = "Threshold")

# Make the labels unique (avoid repetitions) and add the legend
handles, labels = ax[0,0].get_legend_handles_labels() # Get labels
unique = dict(zip(labels, handles)) # Make them unique
ax[0,0].legend(unique.values(), unique.keys(), shadow = 1, facecolor = "whitesmoke") # Add the legend


# q-VALUE SPEARMAN DISTRIBUTION   

# Define an histogram and extract its values
n, bins, rect = ax[0,1].hist(s_q, bins = 30)

# Assign colors based on bin values
for rect, bin_val in zip(rect, bins): # For each rectangle (patch) and value of this rectangle
    if bin_val < threshold: # If the value is lower than the threshold
        rect.set_facecolor('forestgreen') # Significant
        rect.set_label("Significant")
    else: # Otherwise
        rect.set_facecolor('grey')  # No
        rect.set_label("Non-significant")

# Axes settings
ax[0,1].set_title("Spearman q-values Distribution", fontsize = 10)
ax[0,1].set_xlabel("Spearman q-val", fontsize = 10)
ax[0,1].set_ylabel("Frequency", fontsize = 10)
ax[0,1].tick_params(labelsize = 8)
ax[0,1].axvline(0.05, color = 'gray', linewidth = 1, linestyle = '--', label = "Threshold")

# Make the labels unique (avoid repetitions) and add the legend
handles, labels = ax[0,0].get_legend_handles_labels() # Get labels
unique = dict(zip(labels, handles)) # Make them unique
ax[0,1].legend(unique.values(), unique.keys(), shadow = 1, facecolor = "whitesmoke") # Add the legend


# q-VALUES CORRELATION 

# Define a function to assign colors to dots for the scatter plots
def f_color_dots(s,f):
    colors_dots = []
    for i in range(0,len(f)):
        results = (f[i] < 0.05, s[i] < 0.05)
        if True in results: # At least 1 significant
            c = "grey"
            if set(results) == {True}: # Both significant
                c = "green"
        else: 
            c = "black" # None significant
        colors_dots.append(c)
        
    return colors_dots # List with the colors for each q-values pair

# Split data by significance
def split_sig(s_q,   # Spearman q-values
              f_q,   # Fisher q-values
              cols): # Colors 
    d = {"s_q":s_q, "f_q":f_q, "colors":cols} # Dictionary with q-values and colors
    datset = pd.DataFrame(d) # Convert it into a dataframe
    black_sub = datset[datset["colors"] == "black"] # Subset with non-significant q-values
    green_sub = datset[datset["colors"] == "green"] # Subset with significant q-values
    grey_sub = datset[datset["colors"] == "grey"]   # Subset with significant q-values in 1 test
    subsets = (black_sub, green_sub, grey_sub)
    return subsets

subsets = split_sig(s_q, f_q, f_color_dots(s_q,f_q))
sig = ("Non-significant", "Significant in both tests", "Significant in only 1 test")

# Plot the 3 subsets in the same axes 
for i in range(0,len(subsets)): 
    ax[1,0].scatter(list(subsets[i]["f_q"]), # X: Fisher q-values
                  list(subsets[i]["s_q"]), # Y: Spearman q-values
                  s = 11, # Size
                  c = list(subsets[i]["colors"]), # Colors
                  label = sig[i]) # The label of the plot is the significance of the dots

# Axes settings and legend
ax[1,0].set_title("Spearman vs Fisher q-values", fontsize = 10)
ax[1,0].set_xlabel("Fisher q-val", fontsize = 10)
ax[1,0].set_ylabel("Spearman q-val", fontsize = 10)
ax[1,0].tick_params(labelsize = 8) 
ax[1,0].axvline(0.05, color = 'gray', linewidth = 1, linestyle = '--', label = "Threshold")
ax[1,0].axhline(0.05, color = 'gray', linewidth = 1, linestyle = '--')
ax[1,0].legend(shadow = 1, facecolor = "whitesmoke")


# q-VALUES CORRELATION (ZOOMED)

# Define subsets of the q-values 
s_sub = list(df["Spearman_qval"][df["Spearman_qval"] < 0.15][df["Fisher_qval"] < 0.07]) # Spearman
f_sub = list(df["Fisher_qval"][df["Spearman_qval"] < 0.15][df["Fisher_qval"] < 0.07]) # Fisher

# Split data by significance
subsets = split_sig(s_sub,f_sub,f_color_dots(s_sub,f_sub))

# Plot the 3 subsets in the same axes 
for i in range(0,len(subsets)): 
    ax[1,1].scatter(list(subsets[i]["f_q"]),      # X: Fisher q-values
                  list(subsets[i]["s_q"]),        # Y: Spearman q-values
                  s = 20,                         # Size (higher, to simulate zoom effect)
                  c = list(subsets[i]["colors"]), # Colors
                  label = sig[i])                 # The label of the plot is the significance of the dots

# Axes settings and legend
ax[1,1].set_title("Spearman vs Fisher q-values (Zoomed)", fontsize = 10)
ax[1,1].set_xlabel("Fisher q-val", fontsize = 10)
ax[1,1].set_ylabel("Spearman q-val", fontsize = 10)
ax[1,1].tick_params(labelsize = 8)
ax[1,1].axvline(0.05, color = 'gray', linewidth = 1, linestyle = '--', label = "Threshold")
ax[1,1].axhline(0.05, color = 'gray', linewidth = 1, linestyle = '--')
ax[1,1].legend(shadow = 1, facecolor = "whitesmoke")
savefig("qvalues.png")

    
#-#-#-#-#-#-#-#-#-#-#-# DISTRIBUTION of CORRELATIONS #-#-#-#-#-#-#-#-#-#-#-#

# Define data
x = [abs(i[0]) for i in spearman_results] # Absolute value of the correlations, as we interested in the magnitude, not the direction
x2 = [i[0] for i in spearman_results]

# Define the plot (Standard distribution)
fig, ax = plt.subplots(1, 2)
ax[0].hist(x2,
        color = "green",
        alpha = 0.8, bins = 50) # Opacity

# Set title and labels
ax[0].set_title("Distribution of Spearman's correlation values", fontsize = 10)
ax[0].set_xlabel("Spearman correlation", fontsize = 10)
ax[0].set_ylabel("Absolute frequency", fontsize = 10)

# Define the plot (Cumulative distribution)
ax[1].hist(x, 
        cumulative = True, # Cumulative historgram 
        color = "green",
        alpha = 0.8) # Opacity

# Set title and labels
ax[1].set_title("Cumulative Distribution", fontsize = 10)
ax[1].set_xlabel("|Spearman correlation|", fontsize = 10)
ax[1].set_ylabel("Absolute frequency", fontsize = 10)

# Define the asbolute correlation values as a dataframe
df_abs = pd.DataFrame(x)

# Limits of the ranges
values = (0, 0.5, 0.7, 0.8, 0.9, 1)

# Calculate the frequency in each limit
for i in range(0,len(values)-1): 
    n = len(list(df_abs[0][df_abs[0] < (values[i+1])][df_abs[0] >= values[i]]))
    print(f"Between {values[i]} and {values[i+1]} there are {n} values.")

savefig("distrib.png")


#-#-#-#-#-#-#-#-#-#-#-# PLOTTING some ALLELIC (GENOMIC) FREQUENCIES of the SNPs #-#-#-#-#-#-#-#-#-#-#-#

# Define a function to exctract the information we are interested in from a variant
def info_SNP(df,n): 
    
    """
    Function to extract the necessary information from the SNP.
    
    Parameters:
        df = dataframe generated with the VCF data.
        n = index of the SV.
    """   
    
    title_SNP = (df.ID.iloc[n]) # Title = ID
    pos = int(df.POS.iloc[n]) # Position
    chro = df.CHROM.iloc[n] # Chromosome
    f_o = df.Fisher_OddsRatio.iloc[n] # Correlation (Odds ratio) value
    f_q = df.Fisher_qval.iloc[n] # Spearman q-value
    s_c = df.Spearman_Correlation.iloc[n] # Correlation value
    s_q = df.Spearman_qval.iloc[n] # Spearman q-value
    
    return title_SNP, chro, pos, f_o, f_q, s_c, s_q

# Define a function for the creation of the subdataframes and random selecting SNPs    
def subdf_select(spearman_q, # True = significant; False = Non-significant
                  fisher, # True = significant; False = Non-significant
                  n_vars): # Number of selected SNPs
    
    # Spearman's significance
    df_sub = df[df["Spearman_qval"] < 0.4] 
    if spearman_q == False: 
        df_sub = df_sub[df_sub["Spearman_qval"] > 0.1] # So Spearman's q-value ranges 0.1 - 0.4
    else: 
        df_sub = df[df["Spearman_qval"] < 0.05]
    
    # Fisher's significance
    df_sub = df_sub[df_sub["Fisher_qval"] < 0.4]
    if fisher == False: 
        df_sub = df_sub[df_sub["Fisher_qval"] > 0.1]
    else: 
        df_sub = df_sub[df_sub["Fisher_qval"] < 0.05]
        
    # Dataframes with the filtered SNPS
    pos = df_sub[df_sub["Spearman_Correlation"] > 0]
    neg = df_sub[df_sub["Spearman_Correlation"] < 0]

    # Indexes selection
    snps_pos = random.sample(range(len(pos)), k = n_vars)  
    snps_neg = random.sample(range(len(neg)), k = n_vars)
    
    return pos, neg, snps_pos, snps_neg

# Apply the function for all the combinations we want to plot
pos_f005, neg_f005, snps_pos_f005, snps_neg_f005 = subdf_select(True,  True,  6) # Spearman and Fisher significant 
pos_f02,  neg_f02,  snps_pos_f02,  snps_neg_f02  = subdf_select(False, True,  3) # Fisher significant 
pos_005,  neg_005,  snps_pos_005,  snps_neg_005  = subdf_select(True,  False, 6) # Spearman significant 
pos_02,   neg_02,   snps_pos_02,   snps_neg_02   = subdf_select(False, False, 3) # None of two tests significant


# Define some necessary lists
GEBVs = [-4.401016364, -3.184106364, -2.452780833, -2.081470833, 3.170856667, 3.466019167, 3.997353636, 5.002573636] # Genomic breeding values
super_info = [] # Empty list to store information about the variants to then generate a supplementary excel

# Define a function to plot each variant
def plot_snps(subdf,     # Subdataframe
              indexes,   # Indexes with the selected variants
              col,       # Column in the figure (multiplier factor)
              row):      # Row in the figure
    
    # Plot each variant of the df defined in the indexes
    for i in range(0, len(indexes)): 
        snp = subdf.iloc[indexes[i]]
        freqs = [float(a) for a in snp["Frequencies(Genomes)"].split(":")] # Exctract frequencies
        ax[row, i+3*col].scatter(GEBVs, # X-axis = GEBVs values
                                 freqs, # Y-axis = allelic frequencies
                                 color = "black") 
        
        # Apply info_SNP function to exctract important information and append it to super_info list
        info_snp = info_SNP(subdf, indexes[i])
        super_info.append([a for a in info_snp])
        title = info_snp[0]
        if title == ".":
            title  = "No ID"
        ax[row, i+3*col].set_title(title) # Add the title to the plot

# First plot: Fisher significant
fig, ax = plt.subplots(3, 6, # Dimensions (3 rows x 6 columns)
                       figsize = (30, 20),
                       sharex = True) # Common X-axis (GEBVs)
plt.subplots_adjust(bottom = 1.3, right = 1.1, top = 1.8) # Aesthetics
fig.text(0.6, 1.25, 'GEBV', ha = 'center', fontsize = 15) # X-axis label
fig.text(0.09, 1.55, 'Allelic frequency (Genomes)', va = 'center', rotation = 'vertical', fontsize = 15) # Y-axis label
fig.text(0.6, 1.85, 'Allelic frequencies of Spearman (and Fisher) significant SNPs', ha = 'center', fontsize = 30) # Title
plot_snps(pos_f005, snps_pos_f005, 0, 0) # Plot each subdataset
plot_snps(neg_f005, snps_neg_f005, 0, 1)
plot_snps(pos_f02, snps_pos_f02, 0, 2)
plot_snps(neg_f02, snps_neg_f02, 1, 2)
fig.savefig("Plots_Spearman_Fisher.jpg", format = "jpg",  bbox_inches = 'tight') # Save the figure
plt.close(fig)

# Second plot: Fisher non-significant
fig, ax = plt.subplots(3,6,figsize = (30,20),sharex = True)
plt.subplots_adjust(bottom = 1.3, right = 1.1, top = 1.8) # Aesthetics
fig.text(0.6, 1.25, 'GEBV', ha = 'center', fontsize = 15) # X-axis label
fig.text(0.09, 1.55, 'Allelic frequency (Genomes)', va = 'center', rotation = 'vertical', fontsize = 15) # Y-axis label
fig.text(0.6, 1.85, 'Allelic frequencies of Spearman (not Fisher) significant SNPs', ha = 'center', fontsize = 30) # Title
plot_snps(pos_005, snps_pos_005, 0, 0) # Plot each subdataset
plot_snps(neg_005, snps_neg_005, 0, 1)
plot_snps(pos_02, snps_pos_02, 0, 2)
plot_snps(neg_02, snps_neg_02, 1, 2)
fig.savefig("Plots_Spearman_NonFisher.jpg", format = "jpg",  bbox_inches='tight') # Save the figure
plt.close(fig)


# Save the information of the SNPs into an Excel
pd.DataFrame(super_info,columns = ["ID","CHROM","POS","FISHER_ODDS","FISHER_Q","SPEARMAN_CORR","SPEARMAN_Q"]).to_excel("InfoSNPs_Plot.xlsx", index = False)