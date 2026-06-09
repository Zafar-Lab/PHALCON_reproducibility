import pandas as pd
import numpy as np
import os

data = os.path.basename(os.getcwd())
final_vcf_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_nonzero_outputInference.vcf"
final_variant_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_non_zero.avinput.variant_function"
multi_anno_file = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_non_zero.hg19_multianno.txt",sep='\t')
final_vcf = np.loadtxt(final_vcf_file,dtype='object')
final_vcf = pd.DataFrame(final_vcf)
variant_function = pd.read_csv(final_variant_file,sep='\t',header=None)
print(final_vcf[1])
print(variant_function[3])
print(final_vcf[1].tolist())
print(variant_function[3].tolist())
final_variants = []
final_variants_df = variant_function.to_numpy()
for i in final_variants_df:
   if i[6] == '-':
      final_variants.append(i[3]-1)
   else:
      final_variants.append(i[3])
   
if final_vcf.shape[0] == variant_function.shape[0]:
   if [int(i) for i in final_vcf[1].tolist()] == final_variants:
      print("entered")
      final_vcf['variant_function'] = variant_function[0]

final_vcf = final_vcf[~final_vcf['variant_function'].isin(['intronic','intergenic','ncRNA_intronic'])]
final_vcf.drop(columns=['variant_function'],inplace=True)

variant_function = variant_function[~variant_function[0].isin(['intronic','intergenic','ncRNA_intronic'])]

final_vcf = final_vcf.to_numpy()
numcells = final_vcf.shape[1]-9
genotype_arr = np.zeros((final_vcf.shape[0],final_vcf.shape[1]-9+2),dtype='object')
for i in range(final_vcf.shape[0]):
  genotype_arr[i][0] = final_vcf[i][0]
  genotype_arr[i][1] = int(final_vcf[i][1])
  for j in range(final_vcf.shape[1]-9):
    if final_vcf[i][9:][j].split(":")[0] == '0/1':
      genotype_arr[i][j+2] = 1
    elif final_vcf[i][9:][j].split(":")[0] == '0/0':
      genotype_arr[i][j+2] = 0
    else:
      print("Something wrong")

df_genotype = pd.DataFrame(genotype_arr)
df_genotype.drop([0],axis=1,inplace=True)
df_genotype.columns=list(range(numcells+1))
print(df_genotype)



numcells = df_genotype.shape[1]-1
df_genotype = df_genotype.to_numpy()
print(df_genotype)

frames=[]

for variant in range(df_genotype.shape[0]):
   site = df_genotype[variant][0]
   for cell in range(1, numcells+1):
      if df_genotype[variant][cell] == 1:
         multi_anno = multi_anno_file[multi_anno_file['Start']==site]
         multi_anno['Tumor_Sample_Barcode'] = "cell"+str(cell)
         frames.append(multi_anno)

pd.concat(frames).to_csv("maf_file_for_cooccurence_and_mutual_exclusivity.txt",index=False,sep='\t')