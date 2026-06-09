# IMPORTS ###########################################################################################################
import sys
import random
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import betabinom
from sklearn.cluster import SpectralClustering
from sklearn.metrics.pairwise import cosine_similarity, polynomial_kernel, euclidean_distances, nan_euclidean_distances, rbf_kernel
import scipy
from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh,eigs
from scipy.linalg import eigvalsh
from numpy import linalg as LA
from sklearn.neighbors import kneighbors_graph
from sklearn.impute import KNNImputer
import time
from sklearn.metrics import f1_score, recall_score, precision_score,silhouette_score
import copy
import os
from ete3 import Tree
import NNI
import SPR
from sklearn.metrics.cluster import adjusted_rand_score
from scipy import linalg
import re
import scanpy as sc
import leidenalg
import matplotlib.pyplot as plt
sys.setrecursionlimit(100000)
import scipy.cluster.hierarchy as sch

import seaborn as sns 
#####################################################################################################################


import pandas as pd
import numpy as np
import os

data = os.path.basename(os.getcwd())
final_vcf_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_nonzero_outputInference.vcf"
final_variant_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_non_zero.avinput.variant_function"
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

df_genotype.set_index([0],inplace=True)
df_genotype.index.name = None
df_genotype.columns = range(df_genotype.shape[1])

clusters = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_inferred_cluster_labels.txt", header=None, names=["Cluster"])
df_genotype = df_genotype.T
df_genotype["Cluster"] = clusters["Cluster"]
print(df_genotype)

df_clustered = df_genotype.groupby("Cluster").first()

# Perform hierarchical clustering on the 11 clusters to get their order
linkage = sch.linkage(df_clustered, method="ward")
dendro = sch.dendrogram(linkage, no_plot=True)

# Extract the order of clusters from the dendrogram
cluster_order = [df_clustered.index[i] for i in dendro["leaves"]]

# Reorder the full dataset (10,000 cells) based on this cluster order
df_sorted = df_genotype.set_index("Cluster").loc[cluster_order].reset_index()

# Drop the Cluster column for clustermap
df_sorted = df_sorted.drop(columns=["Cluster"])
df_sorted

sns.clustermap(df_sorted, method="ward", cmap="coolwarm", figsize=(12, 6))
plt.show()



'''
argParser = argparse.ArgumentParser(prog='PROG')
argParser.add_argument('-s', '--seed', type=int, default = 179215)
args = argParser.parse_args()
data=os.path.basename(os.getcwd())

outputPrefixName = data + "_variant_calls_heatmap_"



 
def generate_seaborn_colors(n): 
    palette = sns.color_palette("husl", n) 
    rgb_colors = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in palette] 
    return rgb_colors 
 

final_vcf_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_nonzero_outputInference.vcf"
final_variant_file = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+data+"/"+data+"_indels_dbsnp_non_zero.avinput.variant_function"
final_vcf = np.loadtxt(final_vcf_file,dtype='object')
final_vcf = pd.DataFrame(final_vcf)
variant_function = pd.read_csv(final_variant_file,sep='\t',header=None)
print(final_vcf[1])
print(variant_function[3])
print(final_vcf[1].tolist())
print(variant_function[3].tolist())
if final_vcf.shape[0] == variant_function.shape[0]:
   if [int(i) for i in final_vcf[1].tolist()] == variant_function[3].tolist():
      print("entered")
      final_vcf['variant_function'] = variant_function[0]

final_vcf = final_vcf[~final_vcf['variant_function'].isin(['intronic','intergenic','ncRNA_intronic'])]
final_vcf.drop(columns=['variant_function'],inplace=True)

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
print(df_genotype)
df_genotype.drop([0],axis=1,inplace=True)
df_genotype.columns=range(numcells+1)
df_genotype.set_index(df_genotype.columns[0],inplace=True)
df_genotype.index.name = None

df_genotype.columns = range(numcells)

print(df_genotype)
#
#print(df_genotype)
cols=range(1,numcells+1)
data = df_genotype.astype(int)


#print(df_genotype[cols])
print(data)
plt.figure(figsize=(15,9))
sns.clustermap(data)
plt.show()


'''

