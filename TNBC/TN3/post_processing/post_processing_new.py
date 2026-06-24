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

from sklearn.metrics.cluster import adjusted_rand_score
from scipy import linalg
import re
import scanpy as sc
#####################################################################################################################

data = 'TN3'

# first removing clonal dbsnp variants

frames = []
dbsnp_database_file_name = '/home/priya/dbsnp_data/dbsnp_138.b37.vcf'
dbsnp_data = np.loadtxt(dbsnp_database_file_name,dtype='object')
dbsnp_df = pd.DataFrame(dbsnp_data)



index_list = []

vcf_file_name = '/home/priya/Documents/priya/MPT post-processing/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_non_zero_outputInference.vcf'
np_obj = np.loadtxt(vcf_file_name,dtype='object')
vcf_df = pd.DataFrame(np_obj,dtype='object')
print(vcf_df)

genotype_arr = np.zeros((np_obj.shape[0],np_obj.shape[1]-9+4),dtype='object')
df_obtained = pd.DataFrame(genotype_arr,dtype='object')
for i in range(np_obj.shape[0]):
    genotype_arr[i][0] = np_obj[i][0]
    genotype_arr[i][1] = int(np_obj[i][1])
    genotype_arr[i][2] = np_obj[i][3]
    genotype_arr[i][3] = np_obj[i][4]
    for j in range(np_obj.shape[1]-9):
        if np_obj[i][9:][j].split(":")[0] == '0/1':
            genotype_arr[i][j+4] = 1
        elif np_obj[i][9:][j].split(":")[0] == '0/0':
            genotype_arr[i][j+4] = 0
        else:
            print("Something wrong")
numCells = df_obtained.shape[1]-4
cols = [i for i in range(4,numCells+4)]
df_obtained['sum'] = df_obtained[cols].sum(axis=1)
root = df_obtained[df_obtained['sum']== numCells]
root_variants = root[1]
print("Root variants",root_variants)



for i in range(np_obj.shape[0]):
    chrm = vcf_df.iloc[i,0]
    pos = vcf_df.iloc[i,1]
    ref = vcf_df.iloc[i,3] 
    alt = vcf_df.iloc[i,4]

    dbsnp_temp = dbsnp_df[dbsnp_df[1] == pos]
    print(dbsnp_temp)
    for j in range(dbsnp_temp.shape[0]):
        #print("chromosome:",'chr'+str(dbsnp_temp.iloc[j,0]))
        if 'chr'+str(dbsnp_temp.iloc[j,0]) == chrm and dbsnp_temp.iloc[j,3] == ref and dbsnp_temp.iloc[j,4] == alt and pos in root_variants:
            index_list.append(i)
vcf_df.drop(index=index_list, inplace=True)
print(vcf_df)



print(vcf_df)
vcf_np = np.array(vcf_df)
with open('/home/priya/Documents/priya/MPT post-processing/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_dbsnp_non_zero_outputInference.vcf',"w") as vcfFile:
    pos = vcf_np.shape[0]
    numCells = vcf_np.shape[1]-9
    vcfFile.write("##fileformat=VCFv4.1\n")
    vcfFile.write("##source=OurAlgo" + "OurAlgo v" + '0' + "." + '1' + "." + '0' + "\n")
    vcfFile.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    vcfFile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n")
    vcfFile.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n")
    vcfFile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
    vcfFile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
    vcfFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcfFile.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
    vcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

    for i in range(1,numCells+1):
        vcfFile.write("\tcell"+str(i))
    vcfFile.write('\n')

    
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(vcf_np[i][0]+'\t')  # Chr
        vcfFile.write(vcf_np[i][1]+'\t')  # Position
        vcfFile.write(vcf_np[i][2]+'\t')
        vcfFile.write(vcf_np[i][3]+'\t')  # Reference allele
        vcfFile.write(vcf_np[i][4]+'\t')  # Alternate allele
        vcfFile.write(vcf_np[i][5]+'\t')
        vcfFile.write(vcf_np[i][6]+'\t')
        vcfFile.write(vcf_np[i][7]+'\t')
        vcfFile.write(vcf_np[i][8]+'\t')
        for cell in range(numCells):
            vcfFile.write(vcf_np[i][cell+9]+'\t')
        vcfFile.write('\n')

# removing variants from normal bulk files (clonal as well as non clonal)
######################################################(clonal)
mpileup = np.loadtxt(data+'_N.mpileup',dtype='str')
print("mpileup shape :",mpileup.shape)
clonal_list = np.array(mpileup[:,1],dtype='str')
print("Clonal mutations list :",clonal_list)
vcf_obtained = np.loadtxt('/home/priya/Documents/priya/MPT post-processing/'+data+'/Hierarchical_prior_FINAL/'+data+'_orig_hierarchical_dbsnp_non_zero_outputInference.vcf',dtype='str')
print("vcf obtained :\n",vcf_obtained)
vcf_df = pd.DataFrame(vcf_obtained)
clonal_vcf = np.array(vcf_df,dtype='object')


sites_to_be_removed = []
for i in range(clonal_vcf.shape[0]):
  if clonal_vcf[i][1] == mpileup[i][1]:
    site = clonal_vcf[i][1]
    pos = clonal_vcf[i][1]
    ref = clonal_vcf[i][3]
    alt = clonal_vcf[i][4]
    seq = mpileup[i][4]
    ref_count = seq.upper().count(ref.upper())
    alt_count = seq.upper().count(alt.upper())
    if ref_count == 0:
      sites_to_be_removed.append(site)
    else:
      if alt_count/ref_count > 0.1:
        sites_to_be_removed.append(site)
print("Sites being removed :",sites_to_be_removed)
clonal_df = pd.DataFrame(clonal_vcf)
post_process_vcf = np.array(clonal_df[~clonal_df[1].isin(sites_to_be_removed)],dtype='object')


print("Final vcf shape :",post_process_vcf.shape)
tag='/home/priya/Documents/priya/MPT post-processing/'+data+'/Hierarchical_prior_FINAL/'+data
lklhd_df = pd.read_csv(tag+'_orig_hierarchical_final_lklhds.tsv',sep='\t',header=None)
final_df_rc = pd.read_csv(tag+'_orig_hierarchical_final_df.tsv',sep='\t',header=None)
genotype = pd.read_csv(tag+'_orig_hierarchical_inferred_genotypes.tsv',sep=',',header=None)

if lklhd_df.shape[1] == final_df_rc.shape[0]:
  lklhd_df.rename(columns =final_df_rc[1], inplace = True)
  genotype.rename(columns = final_df_rc[1], inplace=True)
  print("Renamed the likelihood matrix columns")

columns = np.array(post_process_vcf[:,1],dtype='int')
column_list = list(columns)
list_indices = []
for i in range(len(columns)-1):
  if abs(columns[i+1]-columns[i])<=10:
    list_indices.append(i)
    list_indices.append(i+1)
remain = [i for i in range(len(list(columns))) if i not in list_indices]
final_site_list = []
for i in remain:
  final_site_list.append(column_list[i])
print("Length of final clonal sites remaning after removing 10 base pair threshold things : ",len(final_site_list))
final_df_rc = final_df_rc[final_df_rc[1].isin(final_site_list)]
final_df_rc.reset_index(inplace=True,drop=True)


numCells = post_process_vcf.shape[1]-9
with open(data+"_post_processed_final.vcf","w") as vcfFile:
    vcfFile.write("##fileformat=VCFv4.1\n")
    vcfFile.write("##source=OurAlgo" + "OurAlgo v" + '0' + "." + '1' + "." + '0' + "\n")
    vcfFile.write("##FILTER=<ID=LowQual,Description=\"Low quality\">\n")
    vcfFile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n")
    vcfFile.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for alt alleles\">\n")
    vcfFile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n")
    vcfFile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
    vcfFile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    vcfFile.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
    vcfFile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")

    for i in range(1,numCells+1):
        vcfFile.write("\tcell"+str(i))
    vcfFile.write('\n')

    pos = post_process_vcf.shape[0]
    numCells = post_process_vcf.shape[1]-9
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(post_process_vcf[i][0]+'\t')  # Chr
        vcfFile.write(post_process_vcf[i][1]+'\t')  # Position
        vcfFile.write(post_process_vcf[i][2]+'\t')
        vcfFile.write(post_process_vcf[i][3]+'\t')  # Reference allele
        vcfFile.write(post_process_vcf[i][4]+'\t')  # Alternate allele
        vcfFile.write(post_process_vcf[i][5]+'\t')
        vcfFile.write(post_process_vcf[i][6]+'\t')
        vcfFile.write(post_process_vcf[i][7]+'\t')
        vcfFile.write(post_process_vcf[i][8]+'\t')
        for cell in range(numCells):
          vcfFile.write(post_process_vcf[i][9+cell]+'\t')
        vcfFile.write('\n')


