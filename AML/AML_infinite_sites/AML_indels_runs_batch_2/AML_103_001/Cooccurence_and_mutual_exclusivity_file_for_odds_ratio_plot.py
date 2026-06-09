import pandas as pd
import numpy as np
import os

data = os.path.basename(os.getcwd())
final_vcf_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_2/"+data+"/"+data+"_indels_dbsnp_nonzero_outputInference.vcf"
final_variant_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_2/"+data+"/"+data+"_indels_dbsnp_non_zero.avinput.variant_function"
final_vcf = np.loadtxt(final_vcf_file,dtype='object')
final_vcf = pd.DataFrame(final_vcf)
variant_function = pd.read_csv(final_variant_file,sep='\t',header=None)
variant_function = variant_function.to_numpy()
print(final_vcf[1])
print(variant_function[3])
print(final_vcf[1].tolist())
print(variant_function[3].tolist())



print(final_vcf)
#final_vcf = final_vcf[~final_vcf['variant_function'].isin(['intronic','intergenic','ncRNA_intronic'])]
#final_vcf.drop(columns=['variant_function'],inplace=True)

#variant_function = variant_function[~variant_function[0].isin(['intronic','intergenic','ncRNA_intronic'])]

final_vcf = final_vcf.to_numpy()
numcells = final_vcf.shape[1]-9
genotype_arr = np.zeros((final_vcf.shape[0],final_vcf.shape[1]-9+3),dtype='object')
for i in range(final_vcf.shape[0]):
  genotype_arr[i][0] = final_vcf[i][0]
  genotype_arr[i][1] = variant_function[i][1] +str("_")+str(final_vcf[i][1])
  genotype_arr[i][2] = variant_function[i][0]
  for j in range(final_vcf.shape[1]-9):
    if final_vcf[i][9:][j].split(":")[0] == '0/1':
      genotype_arr[i][j+3] = 1
    elif final_vcf[i][9:][j].split(":")[0] == '0/0':
      genotype_arr[i][j+3] = 0
    else:
      print("Something wrong")

df_genotype = pd.DataFrame(genotype_arr)
df_genotype = df_genotype[~df_genotype[2].isin(['intronic','intergenic','ncRNA_intronic'])]
df_genotype.drop([0,2],axis=1,inplace=True)
df_genotype.columns=range(df_genotype.shape[1])
#df_genotype.columns=["cell"+str(i) for i in range(numcells+1)]
print(df_genotype)


print(df_genotype.T)
df_genotype.T.to_csv("input_"+data+"_for_odds_ratio_plot.csv",sep=',',header=None,index=False)


