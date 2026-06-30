################ IMPORTS ###################################################################################################################################
import numpy as np
from sklearn.metrics import f1_score, recall_score, precision_score
from sklearn.metrics.cluster import adjusted_rand_score
import pandas as pd
from tqdm import tqdm

data_no = 1
benchmark = 'ADO_decreased_01'

######################################################## TRE Calculation ####################################################################################
def getPairwiseShortestDistance(df_orig,df_obtained):
    numCells = df_orig.shape[0]
    df_orig = df_orig.to_numpy()
    df_obtained = df_obtained.to_numpy()
    cell_pairs = [(i,j) for i in range(numCells) for j in range(numCells) if i<j]
    dist = 0
    for cell_pair in tqdm(cell_pairs,desc="Computing pair wise distances between cells.."):
        cell_1 = cell_pair[0]
        cell_2 = cell_pair[1]
        output_orig = np.add(df_orig[cell_1],df_orig[cell_2])
        output_obtained = np.add(df_obtained[cell_1],df_obtained[cell_2])
        dist+=abs(np.count_nonzero(output_orig==1) - np.count_nonzero(output_obtained==1))
    return dist/(len(cell_pairs))


# Obtained files ##############################################################################################################################################
inferredGenotypeFileName='sc_2000_inferred_genotypes.tsv'   # Inferred genotypes
best_labels = np.loadtxt('/home/priya/Downloads/Final Stage/Benchmarking/'+benchmark+'/Data_'+str(data_no)+'/src/sc_2000_inferred_cluster_labels.txt')   # Inferred cluster labels
df_obtained_rc = pd.read_csv('/home/priya/Downloads/Final Stage/Benchmarking/'+benchmark+'/Data_'+str(data_no)+'/src/sc_2000_final_df.tsv',header=None,sep='\t')     # Final inferred read count matrix


# Ground Truth ##############################################################################################################################################
trueGenofile = '/home/priya/Downloads/Final Stage/Benchmarking/'+benchmark+'/Data_'+str(data_no)+'/simulation/Genotype_data_2000_50.tsv'    # True genotype 
trueLabelfile = '/home/priya/Downloads/Final Stage/Benchmarking/'+benchmark+'/Data_'+str(data_no)+'/simulation/ccm_2000_50.txt'   # True cell to cluster map
trueVCF = '/home/priya/Downloads/Final Stage/Benchmarking/'+benchmark+'/Data_'+str(data_no)+'/simulation/sc_2000_inclusion_output.vcf'



# INFERENCE ##################################################################################################################################################
true_cluster_labels = np.loadtxt(trueLabelfile)  # True cluster labels
true_cluster_labels = np.delete(true_cluster_labels, (0), axis=0)
true_cluster_labels = true_cluster_labels[0]

df_orig = pd.read_csv(trueGenofile,sep='\t',header=None)
True_vcf = pd.read_csv(trueVCF,sep='\t',header=None)
df_orig[0] = True_vcf[1]
df_obtained = pd.read_csv(inferredGenotypeFileName,header=None,index_col=None)
df_obtained = df_obtained.T
df_obtained.insert(loc = 0, column='index',value =df_obtained_rc[1])

mut_pos = list(True_vcf[1])
mut_pos_obtained = list(df_obtained['index'])

for i in mut_pos:
    if i not in mut_pos_obtained:
        l = [i]
        l.extend(list(np.zeros(df_obtained.shape[1]-1)))
        df_obtained.loc[len(df_obtained.index)] = l

rslt_df = df_obtained.sort_values(by = 'index')
rslt_df.drop(columns='index',inplace=True)
df_obtained = rslt_df.T
df_orig.drop(columns=0,inplace=True)
df_orig = df_orig.T

print("Original genotype dataframe shape: ",df_orig.shape)
print("Obtained genotype dataframe shape: ",df_obtained.shape)

print("Computing scores...")
final_ari = adjusted_rand_score(true_cluster_labels,best_labels)
f1 = f1_score(df_orig,df_obtained,pos_label=1,average='micro')
recall = recall_score(df_orig,df_obtained,pos_label=1,average='micro')
precision = precision_score(df_orig,df_obtained,pos_label=1,average='micro')
TRE = getPairwiseShortestDistance(df_orig,df_obtained)

print("Final ARI : ",final_ari)
print("Final f1 score : ",f1)
print("Final recall : ",recall)
print("Final precision : ",precision)
print('Final TRE : ',TRE)