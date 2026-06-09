##################################################################################
# THIS LOOMTOREADCOUNT DOES NOT REMOVE INDELS, TAKES THE MAX ALTERNATE DEPTH OUT OF THE DUPLICATE POSITIONS
# FOR CASES WITH WEIRD ALTERNATE OR REFERENCE COUNTS, IT RANDOMLY SELECTS A REFERENCE AND ALTERNATE ALLELE AND JUST PUTS REFERENCE COUNT AND ALETRNATE COUNTS
# LASTLY, IF THE COVERAGE IS STILL NOT COMPLETE, A COMPLETELY RANDOM BASE IS CHOSEN AND THE REMAINING COVERAGE IS PUT THERE.
##################################################################################
import sys
import random
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import betabinom
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import scipy
from scipy.sparse import csgraph
from scipy.sparse.linalg import eigsh,eigs
from scipy.linalg import eigvalsh
from sklearn.manifold import TSNE
import seaborn as sns
from sklearn.impute import KNNImputer
import time
import loompy
from sklearn.neighbors import kneighbors_graph
from sklearn.metrics import f1_score, recall_score, precision_score
from sklearn.neighbors import NearestCentroid
import copy
import os
import csv
import math
from decimal import *
plt.switch_backend('agg')
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
from sklearn.metrics.cluster import adjusted_rand_score
import re

def indexToChar(index):
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'

def charToIndex(char):
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3

argParser = argparse.ArgumentParser(prog='PROG')

argParser.add_argument('-D', '--dataset', type=str)

args = argParser.parse_args()

data = args.dataset
#loom_name = re.sub('_','-',data)

loomFilename = 'AML-67-001_PDX1.loom'
#mpileupFileName = "AML_84_001.mpileup"
rcFileName = "ReadCounts_"+data+".tsv"
quality_file = "GQ_"+data+".tsv"


with loompy.connect(loomFilename) as ds:
        print('ds layers: ',ds.layers)
        print("ds shape: ",ds.shape)
        gt=ds[''][:]
        print("gt shape: ",gt.shape)
        dp=ds.layers.DP[:]      
        print("dp shape: ",dp.shape)
        gq = ds.layers.GQ[:]
        ad=ds['AD'][:]
        
        print("ad shape: ",ad.shape)
        ro=ds['RO'][:]
        print("ro shape: ",ro.shape)

        rows = ds.shape[0]
        cols = ds.shape[1]
        print("rows: ",rows)
        print("cols: ",cols)
        ad = pd.DataFrame(ad)
        gq = pd.DataFrame(gq)
        dp = pd.DataFrame(dp)
        ro = pd.DataFrame(ro)
        snp=pd.DataFrame({'CHROM':ds.ra['CHROM'],'POS':ds.ra['POS'],'REF':ds.ra['REF'],'ALT':ds.ra['ALT']})
        
        ad_df = pd.DataFrame(ad, index=None)
        ro_df = pd.DataFrame(ro, index=None)
        dp_df = pd.DataFrame(dp, index=None)
        gq_df = pd.DataFrame(gq, index=None)
        gt_df = pd.DataFrame(gt,index=None)

        gq_df = pd.concat([snp, gq_df], axis=1)
        ad_df = pd.concat([snp,ad_df],axis=1)    
        dp_df =pd.concat([snp,dp_df],axis=1)
        ro_df = pd.concat([snp,ro_df],axis=1)
    


print("Finding sites with highest total alternate depth ...")
all_pos_list = list(gt_df.index)

ad_df['sum'] = ad_df.loc[:, 0:dp.shape[1]-1].sum(axis=1)
# group by 'POS' and only keep rows with highest sum

ad_df = ad_df.loc[ad_df.groupby('POS',sort=False)['sum'].idxmax()]
# drop 'sum' column
ad_df = ad_df.drop(columns=['sum'])

kept_indices = ad_df.index


gq_df = gq_df.loc[kept_indices]
dp_df = dp_df.loc[kept_indices]
ro_df = ro_df.loc[kept_indices]
gt_df = gt_df.loc[kept_indices]
snp = snp.loc[kept_indices]

#print("gt df shape :",gt_df.shape)
#print("snp df shape ",snp.shape)


gt_np = gt_df.to_numpy()

gt_df = pd.DataFrame(gt_np,index=None)
gt_df.index = "chr"+snp["CHROM"].map(str) + "_"+snp["POS"].map(str)+"_"+snp["REF"].map(str)+"_"+snp["ALT"].map(str)

#print("snp:",snp)
#print("gt df after keeping max alternate depth instances\n",gt_df)



ad_df = ad_df.drop(columns = ['CHROM','POS','REF','ALT'])
gq_df = gq_df.drop(columns = ['CHROM','POS','REF','ALT'])
ro_df = ro_df.drop(columns = ['CHROM','POS','REF','ALT'])
dp_df = dp_df.drop(columns = ['CHROM','POS','REF','ALT'])

gq = gq_df.to_numpy(dtype='object')
ad = ad_df.to_numpy(dtype='object')
ro = ro_df.to_numpy(dtype='object')
dp = dp_df.to_numpy(dtype='object')

#print(gq)
#print(ad)
#print(ro)
#print(dp)

#print(gq.shape)
#print(ad.shape)
#print(ro.shape)
#print(dp.shape)

#print(gt_df)


all_cols = gt_df.index.to_numpy()

print('Number of sites: ',len(all_cols))

all_rows = gt_df.columns.to_numpy()
print('Number of cells: ',len(all_rows))

#mpileupFile = open(mpileupFileName,'w')

rcFile = open(rcFileName,'w')
gqFile = open(quality_file,'w')
pos_ind=0
c=0
prev = ''
for col in tqdm(all_cols,desc="Constructing readcount file from loom file.."):
    pos_details = col.split('_')
    chr_name = pos_details[0]
    position = pos_details[1]
    
    
    
    ref_allele = pos_details[2]
    alt_allele = pos_details[3]
   

    rcFile.write(chr_name+"\t" + position +"\t" + ref_allele+ "\t"+alt_allele)
    gqFile.write(position)
    for r in range(len(all_rows)):

        alt_depth = ad[c][r]
        coverage = dp[c][r]
        ref_depth = ro[c][r]
        gq_qual = gq[c][r]

        if coverage == 0:
            
            if ref_allele in ['A','C','G','T'] and alt_allele in ['A','C','G','T']:
                rcFile.write("\t"+"0,0,0,0")
            elif ref_allele not in ['A','C','G','T'] or alt_allele not in ['A','C','G','T']:
                rcFile.write("\t"+"0,0")
            else:
                print("Something wrong at coverage = 0 !")
                exit(1)
            gqFile.write("\t"+"0")
        
        else:
            if ref_allele in ['A','C','G','T'] and alt_allele in ['A','C','G','T']:
                counts = np.zeros(4)
                counts[charToIndex(ref_allele)] = int(ref_depth)
                counts[charToIndex(alt_allele)] = int(alt_depth)
            elif ref_allele not in ['A','C','G','T'] or alt_allele not in ['A','C','G','T']:
                counts = np.zeros(2)
                counts[0] = int(ref_depth)
                counts[1] = int(alt_depth)
            else:
                print("Something wrong with indels!")
                exit(1)


            if len(counts)==4 and coverage > ref_depth+alt_depth:
                diff = coverage-(ref_depth+alt_depth)
                other_allele = indexToChar(random.randint(0,3))
                while other_allele == ref_allele or other_allele == alt_allele: 
                    other_allele = indexToChar(random.randint(0,3))
                counts[charToIndex(other_allele)] = int(diff)
            
            rcFile.write('\t')
            if len(counts) == 4:
                rcFile.write(str(int(counts[0]))+','+str(int(counts[1]))+','+str(int(counts[2]))+','+str(int(counts[3])))
            elif len(counts)==2:
                rcFile.write(str(int(counts[0]))+','+str(int(counts[1])))
            else:
                print("Counts wrong!")
                exit(1)
                

            gqFile.write("\t"+str(gq_qual))

    rcFile.write('\n')
    gqFile.write('\n')
    c+=1
    
#mpileupFile.close()
rcFile.close()
gqFile.close()

