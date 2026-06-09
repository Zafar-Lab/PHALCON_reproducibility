import re
import os
import pandas as pd
import numpy as np



folder = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/'
sub_folders = ['AML_12_001']

for aml_file_name in sub_folders:

    os.chdir(folder+aml_file_name+"/")

    current_vcf_file = np.loadtxt(aml_file_name + "_indels_dbsnp_nonzero_weird_variant_outputInference.vcf",dtype='str')

    current_vcf_file_df = pd.DataFrame(current_vcf_file)
    remove_pos_phalcon = [28592635,28592636]
    for pos in remove_pos_phalcon:
        current_vcf_file_df.drop(current_vcf_file_df.loc[current_vcf_file_df[1] == str(pos)].index, inplace=True)
        current_vcf_file_df.reset_index(drop=True, inplace=True)
        
    dbsnp_nz_vcf_np = current_vcf_file_df.to_numpy()
    print("Updated list of variants :",dbsnp_nz_vcf_np)


    with open(aml_file_name + "_indels_dbsnp_nonzero_outputInference.vcf","w") as vcfFile:
            vcfFile.write("##fileformat=VCFv4.1\n")
            vcfFile.write("##source=OurAlgo" + "OurAlgo v" + '0' + "." + '1' + "." + '0' + "\n")
            vcfFile.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
            vcfFile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n")
            vcfFile.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n")
            vcfFile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
            vcfFile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
            vcfFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcfFile.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
            vcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

            pos = dbsnp_nz_vcf_np.shape[0]
            numCells = dbsnp_nz_vcf_np.shape[1]-9

            for i in range(1,numCells+1):
                vcfFile.write("\tcell"+str(i))
            vcfFile.write('\n')

            for i in range(pos):
                vcfFile.write(dbsnp_nz_vcf_np[i][0]+'\t')  # Chr
                vcfFile.write(dbsnp_nz_vcf_np[i][1]+'\t')  # Position
                vcfFile.write(dbsnp_nz_vcf_np[i][2]+'\t')
                vcfFile.write(dbsnp_nz_vcf_np[i][3]+'\t')  # Reference allele
                vcfFile.write(dbsnp_nz_vcf_np[i][4]+'\t')  # Alternate allele
                vcfFile.write(dbsnp_nz_vcf_np[i][5]+'\t')
                vcfFile.write(dbsnp_nz_vcf_np[i][6]+'\t')
                vcfFile.write(dbsnp_nz_vcf_np[i][7]+'\t')
                vcfFile.write(dbsnp_nz_vcf_np[i][8]+'\t')
                for cell in range(numCells):
                    vcfFile.write(dbsnp_nz_vcf_np[i][cell+9]+'\t')
                vcfFile.write('\n')
            
        


