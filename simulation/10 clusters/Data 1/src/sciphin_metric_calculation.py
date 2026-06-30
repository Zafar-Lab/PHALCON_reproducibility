import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, recall_score, precision_score

# output obtained by running sciphin
df_obtained_t = pd.read_csv('/content/drive/MyDrive/Data for comparison/SciPhiN results/result_mut2Sample.tsv',sep='\t')
#dropping the chromosome column
df_obtained_t.drop(['chrom'], axis=1,inplace=True)
# dropping the last column (because of the nans, will figure this out later)
df_obtained_t.drop(df_obtained_t.columns[-1],axis=1,inplace=True)
# sorting the obtained dataframe rows (based on the position)
df_obtained_t.sort_values(by = df_obtained_t.columns[0],inplace=True,ignore_index=True)
df_obtained_t.columns=range(2001)
df_obtained = df_obtained_t.drop_duplicates(subset=0, keep="first")
df_obtained.reset_index(inplace=True)
df_obtained  = df_obtained.drop(['index'],axis=1)
df_for_comparison = np.zeros((df_obtained.shape[0],df_obtained.shape[1]),dtype='int')  # the input df should be in (site x cell) form
df_for_comparison = pd.DataFrame(df_for_comparison)
df_for_comparison[0] = df_obtained[0]
np_for_comparison = np.array(df_for_comparison,dtype = 'int')

# the original vcf matrix, needed for finding out the positions that are mutated
vcf_orig = pd.read_csv('/content/drive/MyDrive/Data for comparison/simulation/sc_2000_inclusion_output.vcf',sep='\t',header=None)
# the original genotype matrix
df_orig = pd.read_csv('/content/drive/MyDrive/Data for comparison/simulation/Genotype_data_2000_50.tsv',sep='\t',header=None)
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
#dropping the first column since the indices are now sorted
df_for_comparison.drop(df_for_comparison.columns[0],axis=1,inplace=True)

df_for_comparison.T.to_csv('Dataframe for comparing results.csv',index=False,header=False)
df_orig = pd.read_csv('Dataframe for comparing results.csv',header=None,index_col=None)
# now df_orig is the extended dataframe containing those sites as well which are their in the obtained tsv file by sciphin. The actually mutated sites are given the original genotypes and the rest are given zero
df_obtained.drop(df_obtained.columns[0],axis=1,inplace=True)
df_obtained.rename(columns={x:y for x,y in zip(df_obtained.columns,range(0,len(df_obtained.columns)))},inplace=True)
df_obtained.T.to_csv('/content/drive/MyDrive/Data for comparison/SciPhiN results/Modified_form_of_obtained_df',index=False,header=False)
df_obtained = pd.read_csv('/content/drive/MyDrive/Data for comparison/SciPhiN results/Modified_form_of_obtained_df',header=None,index_col=None)
f1 = f1_score(df_orig,df_obtained,pos_label=1,average='micro')
recall = recall_score(df_orig,df_obtained,pos_label=1,average='micro')
precision = precision_score(df_orig,df_obtained,pos_label=1,average='micro')
print("F1 score : ",f1)
print("Recall : ",recall)
print("Precision : ",precision)