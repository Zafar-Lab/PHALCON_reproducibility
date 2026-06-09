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
import ast
#####################################################################################################################


argParser = argparse.ArgumentParser(prog='PROG')

argParser.add_argument('-d','--Dataset',type=str)
argParser.add_argument('-i', '--inReadCountFileName', type=str)
argParser.add_argument('-g', '--inQualityFileName', type=str)
argParser.add_argument('-o', '--outputPrefix', type=str)
argParser.add_argument('-r', '--minReadDepth', type=int, default=10)
argParser.add_argument('-q', '--minGenotypeQuality',type=int, default=30)
argParser.add_argument('-a', '--minAltAlleleFrequency', type=float, default=0.2)
argParser.add_argument('-v', '--minVarThreshold', type=float, default=0.5)  
argParser.add_argument('-m', '--minMutantFraction', type=float, default=0.01)
argParser.add_argument('-e', '--eigenValueThreshold', type=float, default=0.6) 
argParser.add_argument('-p', '--maxDecreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-l', '--maxIncreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-t', '--enableGenotypeQualityFilter', type=bool, default=True)
argParser.add_argument('-b', '--treeIterations',type=int,default=30)
argParser.add_argument('-s', '--seed', type=int, default = 120991)
argParser.add_argument('-S','--siteNotPresentInBulk',type=int)


args = argParser.parse_args()
site_not_found_in_bulk = args.siteNotPresentInBulk
data = args.Dataset


def indexToChar(index):
# Converts index to character
    if (index == 0):
        return 'A'
    elif (index == 1):
        return 'C'
    elif (index == 2):
        return 'G'
    elif (index == 3):
        return 'T'



def charToIndex(char):
# Converts character to index
    if (char == 'A'):
        return 0
    elif (char == 'C'):
        return 1
    elif (char == 'G'):
        return 2
    elif (char == 'T'):
        return 3



def charToNum(t):
# Returns int value of an integer string
    if t.isdigit():
        return int(t)
    return t



def getSortedChrPos( keys ):
# Sorting positions for sending it to the VCF file
    key_func = lambda txt: [charToNum(c) for c in re.split('([0-9]+)', txt)]
    return sorted(keys, key = key_func)





def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))



def genoQualityFilter(df,gq,qualityVal):
# Genotype quality filter 
# Puts the readcount as "0,0,0,0" at places where the genotype quality is lower than a threshold
    muts = df.shape[0]
    cells = df.shape[1]-4
    for i in tqdm(range(muts),desc="Applying GENO QUALITY filter..."):
        for j in range(4,4+cells):
            if gq[i][j-4] < qualityVal:
                df[i][j] = "0,0,0,0"



def readDepthFilter(df,minReadDepth):
# Read depth quality filter
# Puts the readcount as "0,0,0,0" at places where the read depth quality is lower than a threshold
    muts = df.shape[0]
    cells = df.shape[1]-4
    for i in tqdm(range(muts),desc="Applying READ DEPTH QUALITY filter..."):
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        read_depth = np.sum(counts,axis=1)
        for j in range(4,4+cells):
            if read_depth[j-4] < minReadDepth:
                df[i][j] = "0,0,0,0"



def altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq):
# Alternate Allele frequency filter (To remove FP due to WGA)
# If there is no alternate allele, put "0,0,0,0" at all cells across the readcount dataframe
# Else, if likelihood > 0.5 and alternate/total depth is less than a certain threshold 
#       then put "0,0,0,0" at that place in the readcount dataframe
    muts = df.shape[0]
    cells = df.shape[1]-4
    pos_retained_list = []

    for i in tqdm(range(muts),desc="Applying ALTERNATE ALLELE FREQUENCY filter..."):
        allele_counts = np.zeros((cells,4),dtype=int)
        ref = df[i][2]
        ref_ind = charToIndex(ref)
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)  # Total allele count of each base at a certain site. Shape : (4,)
        ref_count = allele_counts[ref_ind]
        allele_counts[ref_ind] = -1
        alt_allele = 'X'
        if np.max(allele_counts) > 0:
            alt_ind = np.argmax(allele_counts)
            alt_allele = indexToChar(alt_ind)
            pos_retained_list.append(i)
        if alt_allele=='X':      
            df[i][4:] = np.full(cells, '0,0,0,0')
            continue
        else:
            for cell in range(cells):
                counts = np.array(list((map(int,df[i][4+cell].strip().split(",")))))
                if df[i][4+cell] == '0,0,0,0':
                    continue
                read_depth = sum(counts)
                alt_depth = counts[alt_ind]
                if lklhd_df[i][cell]>0.5 and np.float16(alt_depth/read_depth) < minAltAlleleFreq:
                    df[i][4+cell] = '0,0,0,0'
    return pos_retained_list
                    
                

def variantRemovalFilter(df,threshold):
# Variant removal filter (based on read count information across all cells)
# If read count is  not "0,0,0,0" across more than a certain threshold among all cells then retain that site
    muts = df.shape[0]
    cells = df.shape[1]-4
    pos_retained = []

    for i in tqdm(range(muts),desc="Applying VARIANT REMOVAL filter..."):
        cell_count = 0
        for cell in range(cells):
            if df[i][4+cell]!='0,0,0,0':
                cell_count+=1
        if cell_count >= (threshold*cells):
            pos_retained.append(i)
    print('num of pos retained after third filter: ',len(pos_retained))
    return pos_retained



def variantRemovalMutated(df,lklhd_df,minMutFraction,alt_alleles,pos_retained):
# Variant removal filter (based on fraction of cells mutated)
# If the likelihood is greater than 0.5 and read count is NOT "0,0,0,0" then increase the cell count by 1
#    If cell count is greater than a certain threshold then retain that site
    print('Read Count data dimension at start of VariantRemovalMutated filter: ',df.shape)
    print("GQ dimension at start of VariantRemovalMutated filter: ",gq.shape)
    for i in tqdm(range(df.shape[0]),desc="Applying VARIANT REMOVAL MUTATED filter..."):
        cell_count = 0
        for cell in range(cells):
            if lklhd_df[i][cell]>0.5 and df[i][4+cell]!='0,0,0,0':
                cell_count+=1
        if cell_count >= (minMutFraction*cells):
            pos_retained.append(i)
    print("Length of positions retained after FINAL FILTER: ",len(pos_retained))
    for pos in pos_retained:
        print(df[pos][0],"\t",df[pos][1],"\t",df[pos][2],"\t",' alt inferred: ',alt_alleles[pos],' true_alt: ',df[pos][3])



iterations = args.treeIterations
countFileName = '/home/priya/Documents/priya/AML_indels_batch_2_readcounts_gq_files/'+data+'/ReadCounts_'+data+'.tsv'
qualityFileName ='/home/priya/Documents/priya/AML_indels_batch_2_readcounts_gq_files/'+data+'/GQ_'+data+'.tsv'
outputPrefixName = data+'_orig_bf_test_'
seed=args.seed
eigenValueThreshold = args.eigenValueThreshold
lt = args.maxDecreaseOfClusterCount
gt = args.maxIncreaseOfClusterCount
isQualityFilteringOn = args.enableGenotypeQualityFilter
inferredGenotypeFileName=outputPrefixName+'inferred_genotypes.tsv'
minReadDepth = args.minReadDepth
minAltAlleleFreq = args.minAltAlleleFrequency
varThreshold = args.minVarThreshold
minMutFraction = args.minMutantFraction
qualityValue = args.minGenotypeQuality

if data in ['AML_04_001','AML_07_001','AML_08_001','AML_29_001','AML_41_001','AML_55_001','AML_64_001','AML_75_001']:
    minReadDepth = 5
    qualityValue = 25
    varThreshold = 0.25

random.seed(seed)
np.random.seed(seed)


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1)  


# # # # # # # # # # # # # # # # # # # # # # FINAL ALGORITHM # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
start_time = time.time()
df_for_sites = pd.read_csv("/home/priya/Documents/priya/AML_indels_nonbulk_evident_exonic_sites_dbsnp_nonzero_statistical_BF_test.tsv",sep='\t')
df_for_sites_of_a_sample = df_for_sites[df_for_sites['Sample_id'] == data]
mutated_sites = int(df_for_sites_of_a_sample['PHALCON_inferred_sites'])


df = pd.read_csv(countFileName, sep="\t",header=None)
print(df)
total_sites = df.shape[0]




df = df[df[1].isin([site_not_found_in_bulk])]
df.to_csv("Bulk_evidence_not_found_sites.tsv",sep='\t',header=False)
print(df)
df = df.to_numpy()
print("Dimensions of the input matrix: ",df.shape)


# # # # # # # # # # # # # # # # # # # # # # FILTERS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
filters_start_time = time.time()
muts = df.shape[0]
cells = df.shape[1]-4
alt_alleles = []
retain_indices = []

for i in tqdm(range(muts),desc='Computing Alternate alleles...'):
    allele_counts = np.zeros((cells,4),dtype=int)
    ref = df[i][2]
    ref_ind = charToIndex(ref)
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    allele_counts = np.sum(counts,axis=0)
    ref_count = allele_counts[ref_ind]
    allele_counts[ref_ind] = -1
    if np.max(allele_counts)<=0:
         continue
    retain_indices.append(i)
    alt_ind = np.argmax(allele_counts)
    alt = indexToChar(alt_ind)
    alt_alleles.append(alt)


df = df[retain_indices,:]
df_copy = copy.copy(df)
alt_alleles_copy = copy.copy(alt_alleles)
print('Read count matrix shape after computing alternate alleles is: ',df.shape)



# .................................Quality Filter..............................................................................................
if isQualityFilteringOn: 
    df_gq = pd.read_csv(qualityFileName,sep='\t',header=None)
    df_gq = df_gq[df_gq[0].isin([site_not_found_in_bulk])]
    gq = df_gq.drop(columns=df_gq.columns[0],axis=1).to_numpy()
    print(" Genotype quality matrix :\n",gq)
    print("\n")
    gq = gq[retain_indices,:]
    print("Genotype Quality shape ",gq.shape)
    print("Read count matrix shape ",df.shape)
    genoQualityFilter(df,gq,qualityValue)
print("-------After Genotype quality filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")




# .................................Read depth Filter..............................................................................................
readDepthFilter(df,minReadDepth)
print("-------After Read depth quality filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Calculating likelihood matrix.............................................................................
muts = df.shape[0]
cells = df.shape[1]-4


lklhd_df = pd.DataFrame(index=range(muts), columns = range(cells))
print("Number of Sites: ",muts)
print("Number of Cells: ",cells)
alt_pos = 0
print('Read count matrix shape: ',df.shape,'\nLikelihood matrix shape: ',lklhd_df.shape,"\nLength of computed alternate alleles: ",len(alt_alleles))
print("\n")
lklhd_df = lklhd_df.to_numpy()
print('Computing likelihoods... ')


fp_rate_across_genome = 6.7 * (10**(-5))

for i in tqdm(range(muts),desc='Computing likelihoods... '):
    alt = alt_alleles[alt_pos]
    alt_pos += 1
    altIndex = charToIndex(alt)
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    alt_depth = counts[:,altIndex]

    f_wt = fp_rate_across_genome * total_sites/mutated_sites
    print(f_wt)

    wild_mean_0 = f_wt
    wild_overdispersion_0 = 100
    alpha_wild_0 = wild_mean_0 * wild_overdispersion_0  # genotype 0
    beta_wild_0 = (1 - wild_mean_0) *  wild_overdispersion_0   # genotype 0
    print(alpha_wild_0)
    print(beta_wild_0)


    
    l0 = betabinom.pmf(alt_depth, coverage,alpha_wild_0,beta_wild_0)
    

    mut_mean_2 = f_wt*3
    mut_overdispersion_2 = 100
    alpha_mut_2 = mut_mean_2 *  mut_overdispersion_2   # genotype 2
    beta_mut_2 = (1 - mut_mean_2) * mut_overdispersion_2   # genotype 2
    
    l2 = betabinom.pmf(coverage - alt_depth, coverage,alpha_mut_2, beta_mut_2)

    print(alpha_mut_2)
    print(beta_mut_2)
    
    mut_mean_1 = 0.5 - f_wt
    mut_overdispersion_1 = 2
    mu = 0.2
    alpha_mut_1 = mut_mean_1 * mut_overdispersion_1   # genotype 1
    beta_mut_1 = (1 - mut_mean_1) * mut_overdispersion_1   # genotype 1
    print(alpha_mut_1)
    print(beta_mut_1)

    l1 = ((mu/2)* l0) + ((mu/2)*l2) + ((1-mu)*betabinom.pmf(alt_depth, coverage, alpha_mut_1, beta_mut_1)) 
    
    
    l1 = np.where(l1<0,0,l1)
    l0 = l0/(l0+l1+l2)
    lklhd_df[i] = np.where(coverage==0,0,1-l0)



print("-------After computing likelihood matrix-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")


print("-------After computing likelihood matrix-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Alternate Allele frequency filter.............................................................................
#pos_retained = altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq)
#df = df[pos_retained,:]
#df_copy = df_copy[pos_retained,:]
#alt_alleles = np.array(alt_alleles)
#alt_alleles = alt_alleles[pos_retained]
#lklhd_df = lklhd_df[pos_retained,:]
#gq = gq[pos_retained,:]
print("-------After Alternate Allele Frequency filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Variant Removal filter.............................................................................
#pos_retained_after_third = variantRemovalFilter(df,varThreshold)
#print('Positions retained after Variant removal filter : ',len(pos_retained_after_third))
#df = df[pos_retained_after_third,:]
#df_copy = df_copy[pos_retained_after_third,:]
#alt_alleles = np.array(alt_alleles)
#alt_alleles = alt_alleles[pos_retained_after_third]
#lklhd_df = lklhd_df[pos_retained_after_third,:]
#gq = gq[pos_retained_after_third,:]
print("-------After Variant removal filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Variant Removal mutated filter.............................................................................
#final_pos_retained = []
#variantRemovalMutated(df,lklhd_df,minMutFraction,alt_alleles,final_pos_retained)
end_time = time.time()
#df = df[final_pos_retained,:]
#df_copy = df_copy[final_pos_retained,:]
#alt_alleles = alt_alleles[final_pos_retained]
#gq = gq[final_pos_retained,:]
print("-------After Variant removal mutated filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("-------Filters applied !-------")
print("Time taken for filtering: ",end_time-filters_start_time)
print("\n")



# .................................Calculating likelihood at the remaninig sites.............................................................................
pd.DataFrame(df).to_csv(outputPrefixName+'modified_df_during_filters.tsv',sep='\t',header=False,index=False)
# df = df_copy
df[:,3] = alt_alleles
#lklhd_df = lklhd_df[final_pos_retained,:]
pd.DataFrame(df).to_csv(outputPrefixName+'final_df.tsv',sep='\t',header=False,index=False)
print('Computing lklhd after filters....')
for i in range(df.shape[0]):
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    lklhd_df[i] = np.where(coverage==0,0,lklhd_df[i])

print("-------Final dimensions after all filters are applied-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# # # # # # # # # # # # # # # # # # # # # # GRAPH LAPLACIAN # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
lklhd_df = lklhd_df.T
data_df = pd.DataFrame(lklhd_df)
print(data_df)
data_df =KNNImputer(n_neighbors=8).fit_transform(X=data_df)
pd.DataFrame(data_df).to_csv(outputPrefixName+'final_sciphin_lklhds.tsv',sep='\t',header=False,index=False)
#graph_laplacian = generate_graph_laplacian(df=data_df, nn=20)
print(data_df)

data_df = pd.read_csv(outputPrefixName+'final_sciphin_lklhds.tsv',sep='\t',header=None)
data_df[1] = 1- data_df[0]
data_df[1] = data_df[1].replace(0,0.000001)
data_df[2] = np.log(data_df[1])


data_df.to_csv(outputPrefixName+'final_sciphin_log_lklhds.tsv',sep='\t',header=False,index=False)

print(pd.DataFrame(data_df).isna().all())

print("Total likelihood :",data_df[2].sum())
with open("Log BF test exonic.txt",'a+') as f:
    f.write(str(data_df[2].sum()))
    f.write("\n")

