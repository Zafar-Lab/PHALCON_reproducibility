
import pandas as pd
import numpy as np
from sklearn.metrics import f1_score, recall_score, precision_score
from ete3 import Tree


# The obtained vcf from sciphin results
vcf_obtained = 'snv.vcf'
data_no = 5 # put the data number here

f = open(vcf_obtained)
geno = []
string = ''
for line in f:
  if line.startswith('#'):
    continue
  else:
    ll = line.rstrip().split('\t')
    split = list(ll[-1])
    split[0] = ll[1]
    one_site_geno = split[:-1]
    one_site = [int(i) for i in one_site_geno]
    geno.append(one_site)

# the obtained genotype matrix from sciphin:
#df_obtained_t = pd.read_csv('/content/drive/MyDrive/VCF to Genotype/data_1_10_cluster_mut2Sample.tsv',sep='\t')
# the original vcf matrix, needed for finding out the positions that are mutated

vcf_orig = pd.read_csv('/home/priya/Downloads/Final Stage/Benchmarking/ADO_increased_03/Data_'+str(data_no)+'/simulation/sc_2000_inclusion_output.vcf',sep='\t',header=None)
# the original genotype matrix from simulation
df_orig = pd.read_csv('/home/priya/Downloads/Final Stage/Benchmarking/ADO_increased_03/Data_'+str(data_no)+'/simulation/Genotype_data_2000_50.tsv',sep='\t',header=None)

genotype_arr = np.array(geno)
print("Genotype array created using snv vcf file :\n",genotype_arr)
print("Shape of genotype arr :",genotype_arr.shape)
numcells = genotype_arr.shape[1]-1
df_obtained = pd.DataFrame(genotype_arr)
df_obtained.columns=range(numcells+1)
df_obtained.sort_values(by = df_obtained.columns[0],inplace=True,ignore_index=True)
df_obtained.reset_index(inplace=True)
df_obtained  = df_obtained.drop(['index'],axis=1)
cols = [i for i in range(1,2001)]
df_obtained['sum'] = df_obtained[cols].sum(axis=1)
df_obtained = df_obtained[df_obtained['sum'] != 1]
df_obtained.reset_index(inplace=True)
df_obtained  = df_obtained.drop(['index'],axis=1)
df_obtained  = df_obtained.drop(['sum'],axis=1)
print(df_obtained.shape)
print(df_obtained)

np_obtained = np.array(df_obtained,dtype = 'int')
mutated_pos = set(vcf_orig[1]) #mutated positions in the original genotype matrix
mutated_absent = set(mutated_pos) - set(np_obtained[:,0])
for i in mutated_absent:
  l = list()
  l.append(i)  #appending the index of mutation
  l.extend(list(np.zeros((numcells))))
  l = np.array(l,dtype='int')
  np_obtained = np.r_[np_obtained,[l]]

df_obtained = pd.DataFrame(np_obtained)
df_obtained.sort_values(by = df_obtained.columns[0],inplace=True,ignore_index=True)


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
print(len(set(df_for_comparison[0]) & set(mutated_pos)))
if len(set(df_for_comparison[0]) & set(mutated_pos)) ==50:
  print("All good") # hence all true mutated positions are in sciphin inferred mutations already
else:
  print("TERMINATE")
#dropping the first column since the indices are now sorted

df_for_comparison.drop(df_for_comparison.columns[0],axis=1,inplace=True)
df_for_comparison.T.to_csv('Dataframe_for_comparing_results_more_than_one_cell',index=False,header=False)
df_orig = pd.read_csv('Dataframe_for_comparing_results_more_than_one_cell',header=None,index_col=None)
# now df_orig is the extended dataframe containing those sites as well which are their in the obtained tsv file by sciphin. The actually mutated sites are given the original genotypes and the rest are given zero
# SAVING THE FOLLOWING TO GET THE TREE STRUCTURE AND HENCE CALCULATE MP3 SIMILARITY ON THE SAME
df_obtained.to_csv('Cell_level_information_phylovar_more_than_one.csv',index=False,header=False)
print("shape of what is being saved;",df_obtained.shape)
# THAT"S IT
df_obtained.drop(df_obtained.columns[0],axis=1,inplace=True)
df_obtained.rename(columns={x:y for x,y in zip(df_obtained.columns,range(0,len(df_obtained.columns)))},inplace=True)
df_obtained.T.to_csv('Modified_form_of_obtained_df_more_than_one',index=False,header=False)
df_obtained = pd.read_csv('Modified_form_of_obtained_df_more_than_one',header=None,index_col=None)

f1 = f1_score(df_orig,df_obtained,pos_label=1,average='micro')
recall = recall_score(df_orig,df_obtained,pos_label=1,average='micro')
precision = precision_score(df_orig,df_obtained,pos_label=1,average='micro')
print("Phylovar results:")
print("F1 :",f1)
print("Recall :",recall)
print("Precision :",precision)

#################################################################################################################################################################