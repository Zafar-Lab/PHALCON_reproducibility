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
from sklearn.metrics import f1_score, recall_score, precision_score
import copy
import os
from ete3 import Tree
import NNI
import SPR
from sklearn.metrics.cluster import adjusted_rand_score
from scipy import linalg
import re




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



inferredGenotypeFileName='sc_2000_inferred_genotypes.tsv'


trueGenofile = "/home/priya/Downloads/Final Stage/Benchmarking/7 clusters/Data 2/simulation/Genotype_data_2000_7.tsv"  #CHANGE HERE
df = pd.read_csv(trueGenofile,sep='\t',header=None,index_col=0)
df_orig = df.T

df_obtained = pd.read_csv(inferredGenotypeFileName,header=None,index_col=None)

TRE = getPairwiseShortestDistance(df_orig,df_obtained)

print('final TRE is: ',TRE)
