import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, recall_score, precision_score
# The obtained vcf from sciphin results
data_no = '3'
vcf_obtained = np.loadtxt('data_cpn_05_3_2000_10.vcf',dtype='str')
vcf_df = pd.DataFrame(vcf_obtained)
vcf_df = vcf_df.drop_duplicates(1, keep='first')
vcf_obtained = np.array(vcf_df)
# the obtained genotype matrix from sciphin:
#df_obtained_t = pd.read_csv('/content/drive/MyDrive/VCF to Genotype/data_1_10_cluster_mut2Sample.tsv',sep='\t')
# the original vcf matrix, needed for finding out the positions that are mutated
vcf_orig = pd.read_csv('/home/priya/Downloads/Final_Stage/Benchmarking/Copy_number_05/Data '+data_no+'/simulation/sc_2000_inclusion_output.vcf',sep='\t',header=None)
# the original genotype matrix from simulation
df_orig = pd.read_csv('/home/priya/Downloads/Final_Stage/Benchmarking/Copy_number_05/Data '+data_no+'/simulation/Genotype_data_2000_50.tsv',sep='\t',header=None)
numcells = vcf_obtained.shape[1]-9
genotype_arr = np.zeros((vcf_obtained.shape[0],vcf_obtained.shape[1]-9+2),dtype='object')
for i in range(vcf_obtained.shape[0]):
  genotype_arr[i][0] = vcf_obtained[i][0]
  genotype_arr[i][1] = int(vcf_obtained[i][1])
  for j in range(vcf_obtained.shape[1]-9):
    if vcf_obtained[i][9:][j].split(":")[0] == '0/1':
      genotype_arr[i][j+2] = 1
    elif vcf_obtained[i][9:][j].split(":")[0] == '0/0':
      genotype_arr[i][j+2] = 0
    else:
      print("Something wrong")

mutated_pos = set(vcf_orig[1]) #mutated positions in the original genotype matrix
not_present = [i for i in set(mutated_pos) if i not in set(genotype_arr[:,1])]

for site in not_present:
    l = []
    l.append('chr1')
    l.append(int(site))
    l.extend([int(j) for j in list(np.zeros(genotype_arr.shape[1]-2))])
    l = np.array(l, dtype='object')
    genotype_arr = np.r_[genotype_arr,[l]]

df_obtained = pd.DataFrame(genotype_arr)

df_obtained.drop([0], axis=1,inplace=True)
df_obtained.columns=range(numcells+1)

df_obtained.sort_values(by = df_obtained.columns[0],inplace=True,ignore_index=True)
cols = [i for i in range(1,numcells+1)]
#df_obtained['sum'] = df_obtained[cols].sum(axis=1)
#df_obtained = df_obtained[df_obtained['sum'] != 1]
#df_obtained = df_obtained[df_obtained['sum'] != 0]
df_obtained.reset_index(inplace=True)
df_obtained  = df_obtained.drop(['index'],axis=1)
#df_obtained  = df_obtained.drop(['sum'],axis=1)
df_for_comparison = np.zeros((df_obtained.shape[0],df_obtained.shape[1]),dtype='int')  # the input df should be in (site x cell) form
df_for_comparison = pd.DataFrame(df_for_comparison)
df_for_comparison[0] = df_obtained[0]
np_for_comparison = np.array(df_for_comparison,dtype = 'int')

# changing the first column of the original df to the positions that are mutated
df_orig[0] = vcf_orig[1]
mutated_pos = set(vcf_orig[1]) #mutated positions in the original genotype matrix
for i in range(df_for_comparison.shape[0]):
  l = list()
  l.append(np_for_comparison[i][0])  #appending the index of mutation
  if np_for_comparison[i][0] in mutated_pos:
    l.extend(list(np.array(df_orig[df_orig[0]==np_for_comparison[i][0]])[0][1:]))
    np_for_comparison[i] = l


df_for_comparison = pd.DataFrame(np_for_comparison)
print(df_for_comparison)
print(len(set(df_for_comparison[0]) & set(mutated_pos)))
if len(set(df_for_comparison[0]) & set(mutated_pos)) ==50:
  print("All good") # hence all true mutated positions are in sciphin inferred mutations already
else:
  print("TERMINATE")


#dropping the first column since the indices are now sorted
df_for_comparison.drop(df_for_comparison.columns[0],axis=1,inplace=True)
df_for_comparison.T.to_csv('Dataframe_for_comparing_results',index=False,header=False)
df_orig = pd.read_csv('Dataframe_for_comparing_results',header=None,index_col=None)
print(df_orig.shape)
# now df_orig is the extended dataframe containing those sites as well which are their in the obtained tsv file by sciphin. The actually mutated sites are given the original genotypes and the rest are given zero
df_obtained.drop(df_obtained.columns[0],axis=1,inplace=True)
df_obtained.rename(columns={x:y for x,y in zip(df_obtained.columns,range(0,len(df_obtained.columns)))},inplace=True)
df_obtained.T.to_csv('Modified_form_of_obtained_df',index=False,header=False)
df_obtained = pd.read_csv('Modified_form_of_obtained_df',header=None,index_col=None)
print("Df obtained from sciphi dimension :",df_obtained.shape)
f1 = f1_score(df_orig,df_obtained,pos_label=1,average='micro')
recall = recall_score(df_orig,df_obtained,pos_label=1,average='micro')
precision = precision_score(df_orig,df_obtained,pos_label=1,average='micro')
print("Sciphin results:")
print("F1 :",f1)
print("Recall :",recall)
print("Precision :",precision)