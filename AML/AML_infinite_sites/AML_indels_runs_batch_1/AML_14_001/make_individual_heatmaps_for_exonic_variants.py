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

import seaborn as sns 
#####################################################################################################################


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
sns.clustermap(df_genotype.astype(int),method='single',row_cluster=False)
plt.show()




