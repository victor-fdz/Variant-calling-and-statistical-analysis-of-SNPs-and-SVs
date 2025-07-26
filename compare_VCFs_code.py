import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2
plt.style.use('bmh')
import sys

# Define the files' paths
file01 = sys.argv[1]
file02 = sys.argv[2]

# Define a function to read these files and convert them into tables
def read_file(file):
    
    # Open it and save its lines
    with open(file) as f: 
        lines = f.readlines()
    
    # Iteratively ccess to the lines
    for i in range(0,len(lines)): 
        if not lines[i].startswith("#"): 
            header = lines[0:i-1] # Header
            cols = lines[i-1][:-1].split("\t") # Names of the columns
            rows = [i[:-1].split("\t") for i in lines[i:]] # Rows (variants)
            break

    # Convert it into a table
    df = pd.DataFrame(rows, columns =cols)
    
    return df

df01 = read_file(file01)
df02 = read_file(file02)


def compare_snps(df1,df2):
    
    # Define empty lists to store unique and coincident variants
    unique1 = []
    coincident = []
    
    for i in range(0,len(df1)):
        snp1 = df1.iloc[i] # Define the SNP in df1

        # Elements that must coincide to consider two SNPs as the same:
        must_coin1 = [snp1.POS, snp1.REF, snp1.ALT, snp1.QUAL, snp1.FILTER, snp1.ABA09619, snp1.ABA09622, snp1.ABA09623, snp1.ABA09951, snp1.ABA09952, snp1.ABA09953, snp1.ABA09954, snp1.ABA09955]

        for i in range(0,len(df2)):
            snp2 = df2.iloc[i] # Define the SNP in df2
            must_coin2 = [snp2.POS, snp2.REF, snp2.ALT, snp2.QUAL, snp2.FILTER, snp2.ABA09619, snp2.ABA09622, snp2.ABA09623, snp2.ABA09951, snp2.ABA09952, snp2.ABA09953, snp2.ABA09954, snp2.ABA09955]

            # Compare snps
            same = True
            for i in range(len(must_coin1)):
                if must_coin1[i] != must_coin2[i]: # If some field is different
                    same = False # It is not the same variant
                    break # It is not necessary to compare more fields

            # After comparing snps
            if same == True: # If this "same" object stills True
                coincident.append(snp1) # They are the same variant
                break # It is not necessary to compare it with more variants 

        # After comparing with all the df2
        if same == False: # If this "same" object stills False
            unique1.append(snp1) # The snp i is unique of the df1
            
            
    # Return the lists
    return unique1, coincident

unique1, coincident1 = compare_snps(df01,df02)
unique2, coincident2 = compare_snps(df02,df01)

pdf_filename = "Venn_VCFs.pdf" # File name

with PdfPages(pdf_filename) as pdf:
    fig, ax = plt.subplots()
    # Plot Venn Diagram
    venn2(subsets = (len(unique1),len(unique2),len(coincident1)),
          set_labels = (file01, file02))
    ax.set_title("VCF comparison")
    pdf.savefig(fig)  
    plt.close(fig)