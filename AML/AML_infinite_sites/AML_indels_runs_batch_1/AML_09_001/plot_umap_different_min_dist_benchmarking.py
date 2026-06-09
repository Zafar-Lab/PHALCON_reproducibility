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


import seaborn as sns 
#####################################################################################################################


argParser = argparse.ArgumentParser(prog='PROG')
argParser.add_argument('-D', '--dataset', type=str)
argParser.add_argument("-m","--minDist",type=float)

argParser.add_argument('-o', '--outputPrefix', type=str)


argParser.add_argument('-s', '--seed', type=int, default = 1719215)

args = argParser.parse_args()


data=os.path.basename(os.getcwd())

min_dist_defined = args.minDist

outputPrefixName = data + "_indels_benchmarking_"




inferredClusterLabelFileName = outputPrefixName+"inferred_cluster_labels"+".txt"

 
def generate_seaborn_colors(n): 
    palette = sns.color_palette("husl", n) 
    rgb_colors = [(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in palette] 
    return rgb_colors 
 
# Example usage 


data_df =pd.read_csv(outputPrefixName+'final_lklhds.tsv',sep='\t',header=None)
cluster = np.loadtxt(inferredClusterLabelFileName)
cluster= cluster.astype(int)
adata = sc.AnnData(X=data_df)
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X')
sc.tl.umap(adata,   min_dist = min_dist_defined)

n= len(set(cluster))
distinct_colors = generate_seaborn_colors(n) 



random_colors=dict()
for i in range(len(set(cluster))):
   random_colors[i] = distinct_colors[i]
colors=dict()
for i in range(len(cluster)):
   colors[i] = random_colors[cluster[i]]
labels = pd.DataFrame(colors.items(), columns=['index', 'label'])
adata.obs['batch'] = cluster.astype('str')
adata.uns['batch_colors'] = random_colors.values()
sc.pl.umap(adata,color='batch',show=False)
plt.savefig(outputPrefixName+"_umap.svg")
plt.savefig(outputPrefixName+"_umap.pdf")


