import pandas as pd
import numpy as np


post_process_vcf = np.loadtxt("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/AML_09_001/AML_09_001_indels_dbsnp_nonzero_outputInference.vcf",dtype='object')
post_process_vcf[4][7] = "SVTYPE=DUP;END=28608317"
# the name "post_process_vcf" can be misleading, i just took this code from somewhere where the code was written

numCells = post_process_vcf.shape[1]-9
with open("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/AML_09_001/AML_09_001_vep_indels_dbsnp_nonzero_outputInference.vcf","w") as vcfFile:
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

    for i in range(1,numCells+1):
        vcfFile.write("\tcell"+str(i))
    vcfFile.write('\n')

    pos = post_process_vcf.shape[0]
    numCells = post_process_vcf.shape[1]-9
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(post_process_vcf[i][0]+'\t')  # Chr
        vcfFile.write(post_process_vcf[i][1]+'\t')  # Position
        vcfFile.write(post_process_vcf[i][2]+'\t')
        vcfFile.write(post_process_vcf[i][3]+'\t')  # Reference allele
        vcfFile.write(post_process_vcf[i][4]+'\t')  # Alternate allele
        vcfFile.write(post_process_vcf[i][5]+'\t')
        vcfFile.write(post_process_vcf[i][6]+'\t')
        vcfFile.write(post_process_vcf[i][7]+'\t')
        vcfFile.write(post_process_vcf[i][8]+'\t')
        for cell in range(numCells):
          vcfFile.write(post_process_vcf[i][9+cell]+'\t')
        vcfFile.write('\n')