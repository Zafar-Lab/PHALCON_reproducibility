import re
import os
import pandas as pd
import numpy as np



folder = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/'
sub_folders = ['AML_38_003']

for aml_file_name in sub_folders:
    print("File working on: ",aml_file_name)

    os.chdir(folder+aml_file_name+"/")
    annovar_file = pd.read_csv(aml_file_name + "_indels_dbsnp_non_zero.avinput.variant_function",sep='\t',header=None)
    annovar_file_np = annovar_file.to_numpy(dtype='object')
    dbsnp_nz_vcf_file_np = np.loadtxt(aml_file_name + "_indels_dbsnp_nonzero_outputInference.vcf",dtype='object',skiprows=10)
  
    dbsnp_nz_vcf_file = pd.DataFrame(dbsnp_nz_vcf_file_np,dtype='object')
    phalcon_variants = dbsnp_nz_vcf_file[1].tolist()


   

    remove_pos_annovar = [148504717]

    remove_pos_phalcon = [148504717,148504716]


    print(remove_pos_annovar)
    print(remove_pos_phalcon)
    if remove_pos_annovar and remove_pos_phalcon:

        print("vcf and annovar file changed in:",aml_file_name)
        #os.system("mv "+aml_file_name + "_indels_dbsnp_non_zero.avinput.variant_function "+aml_file_name + "_indels_dbsnp_non_zero_redundant.avinput.variant_function")
        #os.system("mv "+aml_file_name + "_indels_dbsnp_nonzero_outputInference.vcf "+aml_file_name + "_indels_dbsnp_nonzero_redundant_outputInference.vcf")
   

        for pos in remove_pos_annovar:
            annovar_file.drop(annovar_file.loc[annovar_file[3] == pos].index, inplace=True)
            annovar_file.reset_index(drop=True, inplace=True)
            
        for pos in remove_pos_phalcon:
            dbsnp_nz_vcf_file.drop(dbsnp_nz_vcf_file.loc[dbsnp_nz_vcf_file[1] == str(pos)].index, inplace=True)
            dbsnp_nz_vcf_file.reset_index(drop=True, inplace=True)
        


        annovar_file.to_csv(aml_file_name + "_indels_dbsnp_non_zero.avinput.variant_function",sep='\t',header=False,index=False)
        dbsnp_nz_vcf_np = dbsnp_nz_vcf_file.to_numpy()

        print("Shapes for annotated file and final variants read count are same : ",dbsnp_nz_vcf_np.shape[0]==annovar_file.shape[0])

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
            
        


