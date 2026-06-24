import numpy as np
import pandas as pd
import os


data = 'TN2'


vcf_file_name = '/home/priya/Downloads/Final_Stage/Real Data Results/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_outputInference.vcf'
np_obj = np.loadtxt(vcf_file_name,dtype='object')
vcf_df = pd.DataFrame(np_obj,dtype='object')
print(vcf_df)
genotype_arr = np.zeros((np_obj.shape[0],np_obj.shape[1]-9+4),dtype='object')
df_obtained = pd.DataFrame(genotype_arr,dtype='object')
for i in range(np_obj.shape[0]):
    genotype_arr[i][0] = np_obj[i][0]
    genotype_arr[i][1] = int(np_obj[i][1])
    genotype_arr[i][2] = np_obj[i][3]
    genotype_arr[i][3] = np_obj[i][4]
    for j in range(np_obj.shape[1]-9):
        if np_obj[i][9:][j].split(":")[0] == '0/1':
            genotype_arr[i][j+4] = 1
        elif np_obj[i][9:][j].split(":")[0] == '0/0':
            genotype_arr[i][j+4] = 0
        else:
            print("Something wrong")
numCells = df_obtained.shape[1]-4
cols = [i for i in range(4,numCells+4)]
df_obtained['sum'] = df_obtained[cols].sum(axis=1)
zero = df_obtained[df_obtained['sum']== 0]

print(zero.index)
vcf_df.drop(index=zero.index.tolist(),inplace=True)
print(vcf_df)
final_df_manual = pd.read_csv('/home/priya/Downloads/Final_Stage/Real Data Results/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_final_df_check.tsv',header=None,sep=',')
indices_to_keep = list(final_df_manual[1])
indices_to_keep = [str(i) for i in indices_to_keep]
print(indices_to_keep)
print("vcf df\n",vcf_df)
vcf_df = vcf_df[vcf_df[1].isin(indices_to_keep)]
#zero = df_post.drop(columns = ['sum'], inplace=True)
print(vcf_df)
  #merged_df = vcf_df.merge(zero, on=['CHR','POS', 'REF','ALT'], how='outer', indicator=True)
#print(merged_df)
#result_df = merged_df.drop(columns=['_merge'])
#result_df.drop_duplicates(subset=None,keep=False,inplace=True)

#result_df.rename(columns={'CHR_x':"CHR"},inplace=True)
#result_df.drop(columns=['CHR_y'],inplace=True)

print(vcf_df)
vcf_np = np.array(vcf_df)
with open('/home/priya/Downloads/Final_Stage/Real Data Results/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_non_zero_outputInference.vcf',"w") as vcfFile:
    pos = vcf_np.shape[0]
    numCells = vcf_np.shape[1]-9
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

    
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(vcf_np[i][0]+'\t')  # Chr
        vcfFile.write(vcf_np[i][1]+'\t')  # Position
        vcfFile.write(vcf_np[i][2]+'\t')
        vcfFile.write(vcf_np[i][3]+'\t')  # Reference allele
        vcfFile.write(vcf_np[i][4]+'\t')  # Alternate allele
        vcfFile.write(vcf_np[i][5]+'\t')
        vcfFile.write(vcf_np[i][6]+'\t')
        vcfFile.write(vcf_np[i][7]+'\t')
        vcfFile.write(vcf_np[i][8]+'\t')
        for cell in range(numCells):
            vcfFile.write(vcf_np[i][cell+9]+'\t')
        vcfFile.write('\n')

vcf_file_name = '/home/priya/Downloads/Final_Stage/Real Data Results/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_non_zero_outputInference.vcf'
np_obj = np.loadtxt(vcf_file_name,dtype='object')
vcf_df = pd.DataFrame(np_obj,dtype='object')
print(vcf_df)
genotype_arr = np.zeros((np_obj.shape[0],np_obj.shape[1]-9+4),dtype='object')
df_obtained = pd.DataFrame(genotype_arr,dtype='object')
for i in range(np_obj.shape[0]):
    genotype_arr[i][0] = np_obj[i][0]
    genotype_arr[i][1] = int(np_obj[i][1])
    genotype_arr[i][2] = np_obj[i][3]
    genotype_arr[i][3] = np_obj[i][4]
    for j in range(np_obj.shape[1]-9):
        if np_obj[i][9:][j].split(":")[0] == '0/1':
            genotype_arr[i][j+4] = 1
        elif np_obj[i][9:][j].split(":")[0] == '0/0':
            genotype_arr[i][j+4] = 0
        else:
            print("Something wrong")
numCells = df_obtained.shape[1]-4


gene_name_file = '/home/priya/Downloads/Final_Stage/annotation by annovar/Bam post processed files FINAL/Final Post processed MPT/'+data+'/'+data+'_final.avinput.variant_function'
gene_name = np.loadtxt(gene_name_file,dtype='object')

df_obtained.insert(loc=2,column='gene_name',value='gene')

for i in range(df_obtained.shape[0]):
    chrm = df_obtained.iloc[i,0]
    pos = df_obtained.iloc[i,1]
    gene = gene_name[i][1]
    df_obtained.iloc[i,2] = gene

print("df_obtained final:\n",df_obtained)
cols = [i for i in range(4,numCells+4)]

df_obtained['sum'] = df_obtained[cols].sum(axis=1)
df_obtained.drop(columns = cols, inplace=True)
df_obtained.to_csv('fraction_of_cells_mutated.tsv',sep='\t')
