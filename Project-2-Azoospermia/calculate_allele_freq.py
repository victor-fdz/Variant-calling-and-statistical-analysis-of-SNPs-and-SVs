### This script takes a VCF file and generates a BED file with only the information required + read frequences and putative genotype.
### To execute this script you must install the module "pysam" in your conda environment. 

# Import sys module. First argument in the terminal is the input and the second one is the name of the output.
import sys
vcffile = sys.argv[1]
outfile=sys.argv[2]

# Import pysam module.
import pysam

# Define the function to calculate the allele frequences and the putative genotype.
def allele_freqs_and_genotype(vcf_file):
    
    # Define the VCF file.
    vcf = pysam.VariantFile(vcf_file)
    
    # To append the information of each variant.
    information=[]
    
    # Iterate over each record (variant) in the VCF file.
    for record in vcf:
        
        # Filter out lines that do not have "PASS" in the FILTER column.
        if 'PASS' not in record.filter:
            continue
            
        # Search for the number of reads supporting SV or reference allele.
        rr = record.samples[0].get('RR', 0)
        rv = record.samples[0].get('RV', 0)
        dr = record.samples[0].get('DR', 0)
        dv = record.samples[0].get('DV', 0)
        
        # Search for other important info you could need.
        variantID = record.id
        chromosome = record.chrom
        alternative = record.alts
        start = record.pos
        refe = record.ref
        end = record.stop
        other_chr = record.info.get('CHR2','Not_BND')
        other_pos = record.info.get('POS2','Not_BND')
        gt = record.samples[0].get('GT', 'UNKNOWN')
        genotype = str(gt[0])+'/'+str(gt[1])
        SVtype = record.info.get('SVTYPE','UKNOWN')
        length = end - start
        if SVtype == 'INS':
          length = record.info['INSLEN']
        
        
        # Calculate frequencies.
        preciseness=record.info.get('PRECISE','IMPRECISE') 
        
        # If the variant is PRECISE use junction reads.
        if preciseness:
          preciseness='PRECISE'
          total_reads = rr + rv
          if total_reads > 0:
            variant_frequence = rv / total_reads
            ref_frequence = 1 - variant_frequence
          else:
            ref_frequence = variant_frequence = 0.0
            
        # If is not PRECISE (it is IMPRECISE) use paired-end reads.
        else:
          preciseness='IMPRECISE'
          total_reads = dr + dv
          if total_reads > 0:
            variant_frequence = dv / total_reads
            ref_frequence = 1 - variant_frequence
          else:
            ref_frequence = variant_frequence = 0.0
        
        # If the depth is so low, do not consider the SV.
        if total_reads<8:
          continue
        
        # Putative genotypes depending on the frequences.
        if  1 >= variant_frequence >= 0.9:
          put_genotype='hzALT'
        elif 0.9 > variant_frequence >= 0.7:
          put_genotype='hzALT?'
        elif 0.7 > variant_frequence >= 0.6:
          put_genotype='Het?'
        elif 0.6 > variant_frequence > 0.4:
          put_genotype='Het'
        elif 0.4 > variant_frequence >= 0.3:
          put_genotype='Het?'
        elif 0.3 > variant_frequence > 0.1:
          put_genotype='hzWT?'
        elif 0.1 >= variant_frequence >= 0:
          put_genotype='hzWT'
        
        #Just in case there was an error with the allelic frequences calculation.
        else:
          continue 
            
        # Create a dictionary with the necessary information of the SV and append it to the "information" list.
        info={
        'chr':chromosome,
        'variantID':variantID,
        'precise':preciseness,
        'start':start,
        'stop':end,
        'other_chr':other_chr,
        'other_pos':other_pos,
        'length':length,
        'alternative':alternative,
        'refe':refe,
        'variantType':SVtype,
        'delly_gen':genotype,
        'putative_genotype':put_genotype,
        'ref': ref_frequence,
        'variant': variant_frequence,
        }
        
        information.append(info)

    return information

# Appply the function and write the results into an output file. 
info = allele_freqs_and_genotype(vcffile)
with open (outfile,'w') as outfile:
  outfile.write("#CHR\tSTART\tEND\tCHR2\tPOS2\tLEN\tREFERENCE\tALT\tSV_TYPE\tID\tPRECISE\tDELLY_GENOTYPE\tREF_FREQ\tSV_FREQ\tPUT_GENOTYPE\n")
  for element in info:
    outfile.write(f"{element['chr']}\t{element['start']}\t{element['stop']}\t{element['other_chr']}\t{element['other_pos']}\t{element['length']}\t{element['refe']}\t{element['alternative']}\t{element['variantType']}\t{element['variantID']}\t{element['precise']}\t{element['delly_gen']}\t{element['ref']}\t{element['variant']}\t{element['putative_genotype']}\n")
