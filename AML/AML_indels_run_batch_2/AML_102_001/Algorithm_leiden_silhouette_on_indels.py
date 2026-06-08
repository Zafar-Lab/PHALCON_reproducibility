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
argParser.add_argument('-i', '--inReadCountFileName', type=str)
argParser.add_argument('-g', '--inQualityFileName', type=str)
argParser.add_argument('-o', '--outputPrefix', type=str)
argParser.add_argument('-r', '--minReadDepth', type=int, default=10)
argParser.add_argument('-q', '--minGenotypeQuality',type=int, default=30)
argParser.add_argument('-a', '--minAltAlleleFrequency', type=float, default=0.2)
argParser.add_argument('-v', '--minVarThreshold', type=float, default=0.5)  
argParser.add_argument('-m', '--minMutantFraction', type=float, default=0.01)
argParser.add_argument('-e', '--eigenValueThreshold', type=float, default=0.5) 
argParser.add_argument('-p', '--maxDecreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-l', '--maxIncreaseOfClusterCount', type=int, default=0)
argParser.add_argument('-t', '--enableGenotypeQualityFilter', type=bool, default=True)
argParser.add_argument('-b', '--treeIterations',type=int,default=30)
argParser.add_argument('-s', '--seed', type=int, default = 17215)
argParser.add_argument("-min","--minDist",type=float,default=0.2)

args = argParser.parse_args()



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



def getRightOrderForVCF(df):
# Returns the final list of indices to be sent to the vcf file
    final_indices = []
    unique_chr_keys = np.unique(df[:,0])
    chr_keys = np.array(getSortedChrPos(unique_chr_keys))
    dict={}
    for key in chr_keys:
        dict[key] = []
    for row in range(df.shape[0]):
        dict[df[row][0]].append(row)
    for key in dict.keys():
        for ind in dict[key]:
            final_indices.append(ind)
    return final_indices



def writeVCFHeader(vcfFile):
# Writes the header of the VCF file
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



def writeEntry(df,gq,vcfFile,inferred_genotypes):  #df = Read count dataframe (sites x cells)
# Writes entry of the final positions in the VCF file
    pos = df.shape[0]
    numCells = df.shape[1]-4
    alt_depth = 0
    count_entry = []
    for i in range(pos):
        vcfFile.write(str(df[i][0])+'\t')  # Chr
        vcfFile.write(str(df[i][1])+'\t')  # Position
        vcfFile.write('*\t')
        vcfFile.write(str(df[i][2])+'\t')  # Reference allele
        vcfFile.write(str(df[i][3])+'\t')  # Alternate allele
        vcfFile.write('*\t')
        vcfFile.write('PASS\t')
        vcfFile.write('DP=')
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        depthAllCells = np.sum(allele_counts)
        vcfFile.write(str(depthAllCells)+'\t')
        vcfFile.write('GT:AD:DP:GQ:PL\t')
        for cell in range(numCells):
            if int(inferred_genotypes[i][cell])==1:
                vcfFile.write('0/1:')
                qual = str(gq[i][cell])
            else:
                vcfFile.write('0/0:')
                qual = str(gq[i][cell])
            count_entry = list(map(int,df[i][4+cell].strip().split(",")))
            if len(count_entry)==4:
               alt_depth = int(count_entry[charToIndex(df[i][3])])  # we first find the alternate allele (df[i][3]) adn then we convert it to index and then pick the entry at that particular index
            elif len(count_entry)==2:
               alt_depth = int(count_entry[1])   # indels are different. indels only have ref:alt, i.e. its count will look like [48,19], so we just choose the number written at alt i.e. count_entry[1]
            else:
               print("Error occured while making vcf file")
            vcfFile.write(str(alt_depth)+':')
            total_depth = str(np.sum(np.array(count_entry)))  # sum is the total depth in all cases
            vcfFile.write(total_depth+':'+qual+'\t')
        vcfFile.write('\n')



def writeVCFFile(df,gq,vcfFile,inferred_genotypes):
# Writes the whole VCF file (header + final positions)
	writeVCFHeader(vcfFile)
	writeEntry(df,gq,vcfFile,inferred_genotypes)



def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))



def genoQualityFilter(df,gq,qualityVal):
# Genotype quality filter 
# Puts the readcount as "0,0,0,0" at places where the genotype quality is lower than a threshold
    muts = df.shape[0]
    cells = df.shape[1]-4
    cmm1=0
    #print("gq filter gq \n",gq)
    #print("##### SEE THIS \n at geno filter, df looks like ",df)
    #print("##### SEE THIS \n at geno filter, df[4][2] looks like ",df[4][2])
    #print("##### SEE THIS \n at geno filter, df[4][3] looks like ",df[4][3])
    for i in tqdm(range(muts),desc="Applying GENO QUALITY filter..."):
        for j in range(4,4+cells):
            if gq[i][j-4] < qualityVal:
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   df[i][j] = "0,0,0,0"
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   cmm1+=1
                   df[i][j] = "0,0"
                else:
                   print("Something wrong in gq filter at pos",i," cell ",j)
    print("how many were made to be zero at geno qual filter ",cmm1)



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
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   df[i][j] = "0,0,0,0"
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   df[i][j] = "0,0"
                else:
                   print("Something wrong in read depth filter at pos",i," cell ",j)
  


def altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq):
# Alternate Allele frequency filter (To remove FP due to WGA)
# If there is no alternate allele, put "0,0,0,0" at all cells across the readcount dataframe
# Else, if likelihood > 0.5 and alternate/total depth is less than a certain threshold 
#       then put "0,0,0,0" at that place in the readcount dataframe
    muts = df.shape[0]
    cells = df.shape[1]-4
    pos_retained_list = []

    for i in tqdm(range(muts),desc="Applying ALTERNATE ALLELE FREQUENCY filter..."):
        if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
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
        elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
            allele_counts = np.zeros((cells,2),dtype=int)
            counts = getAltAlleleCounts(df[i][4:])
            counts = np.array(counts.tolist())

            allele_counts = np.sum(counts,axis=0)  # Total allele count of each base at a certain site. Shape : (2,)
          
            alt_allele = 'X'
            allele_counts[0] = -1  # in case of indels, reference index is always zero and alternate index is always 1 because its counts go like ref:alt
            if np.max(allele_counts) > 0:
                alt_ind = np.argmax(allele_counts)
                alt_allele = df[i][3]
                pos_retained_list.append(i)
        else:
            print("Something wrong in alt allele freq filter at pos ",i)

        if alt_allele=='X':      
            if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
               df[i][4:] = np.full(cells, '0,0,0,0')
               continue
            elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
               df[i][4:] = np.full(cells, '0,0')
               continue
            else:
               print("Something wrong with alternate allele at alt allele freq filter at pos ",i)
        else:
            for cell in range(cells):
                counts = np.array(list((map(int,df[i][4+cell].strip().split(",")))))
                if df[i][4+cell] == '0,0,0,0' or df[i][4+cell] == '0,0':
                    continue
                read_depth = sum(counts)
                if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                   alt_depth = counts[alt_ind]
                elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                   alt_depth = counts[1]  # because indels carry readcount info as ref:alt, so reference depth is present at index 0 and alternate depth is present at index 1   

                if lklhd_df[i][cell]>0.1 and np.float16(alt_depth/read_depth) < minAltAlleleFreq:
                    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
                       df[i][4+cell] = '0,0,0,0'
                    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
                       df[i][4+cell] = "0,0"
                    else:
                       print("Something wrong with read count information at later part of alt allele freq filter at pos ",i," and cell ",cell)
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
            if df[i][4+cell] not in ['0,0,0,0','0,0']:
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
            if lklhd_df[i][cell]>0.5 and df[i][4+cell] not in ['0,0,0,0','0,0']:
                cell_count+=1
        if cell_count >= (minMutFraction*cells):
            pos_retained.append(i)
    print("Length of positions retained after FINAL FILTER: ",len(pos_retained))
    for pos in pos_retained:
        print(df[pos][0],"\t",df[pos][1],"\t",df[pos][2],"\t",' alt inferred: ',alt_alleles[pos],' true_alt: ',df[pos][3])



def chunk_ML(splits, inverted_splits, ones, zeros,split_names):
# Calculate the maximum likelihood for each site
    tmp = np.matmul(splits,ones)+np.matmul(inverted_splits,zeros)
    pd.DataFrame(tmp).to_csv("tmp_matrix.tsv",sep='\t',index=False,header=False)
    pd.DataFrame(split_names).to_csv("genotype_name_tmp.tsv",sep='\t',index=False,header=False)
    par_case = 0
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Par':
          par_case+=1

    loss_case = 0
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Loss':
          loss_case+=1

    priors = np.zeros((splits.shape[0],ones.shape[1]))
    for i in range(split_names.shape[0]):
       if split_names[i][-1] == 'Par':
          priors[i] = np.log(  np.exp(-700) * 0.001 * (1/par_case)  ) 
       elif split_names[i][-1] == 'Loss':
          priors[i] = np.log(  np.exp(-700) * 0.001 * (1/loss_case) ) 


    tmp = tmp + priors
    split_indexes = np.argmax(tmp,axis=0)
    genotypes = splits[split_indexes]
    split_names = split_names[split_indexes]
    L = np.sum(np.max(tmp, axis=0))
    return (L, genotypes, split_names)



def heterozygous(ccm, tree, rc, names, split_list, split_counts, split_names,clusters):
   for node in tree.traverse("preorder"):
        if node.is_leaf():
          mutated_names = set([node.name])
        else:
          mutated_names = set(node.get_leaf_names(is_leaf_fn=None))
        tmp_arr = np.zeros(rc.shape[1])
        for name in names:
           if name in list(mutated_names):
              for cell in ccm[name]:
                tmp_arr[cell] = 1
        tmp_arr = tmp_arr.tolist()
        tmp_copy = tmp_arr.copy()
        if tmp_arr not in split_list:
            split_list.append(tmp_arr)
            split_counts.append(tmp_arr.count(1))
            if tmp_arr.count(1) == len(tmp_arr):
               tmp_copy.append('Clonal')
            else:
               tmp_copy.append('Het')
            split_names.append(tmp_copy)
   return split_counts, split_list, split_names

def parallel(ccm, tree, rc, names, split_list, split_counts,split_names,clusters):
   m = []  #list to hold mutated names (includes duplicates)
   for node1 in tree.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  # descendants of current node
     ans = list( h for h in node1.get_ancestors() )  # anscestors of current node
     for node2 in tree.traverse("preorder"):
      if len(ans)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 not in des and node1.up!=node2.up:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             mutated_names2 = set([node2.name])
           else:
             mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m.append( mutated_names1.union( mutated_names2))
   mutated_names = []
   for i in m:
     if i not in mutated_names:
       mutated_names.append(i)    # List of all possibilities of parallel mutation (no duplicate cases)
   for j in mutated_names:
     not_mutated_names = clusters - j
     tmp_arr = np.zeros(rc.shape[1])
     for name in names:
       if name in list(j):
         for cell in ccm[name]:
           tmp_arr[cell] = 1
     tmp_arr = tmp_arr.tolist()
     tmp_copy = tmp_arr.copy()
     if tmp_arr not in split_list:
       split_list.append(tmp_arr)
       split_counts.append(tmp_arr.count(1))
       tmp_copy.append('Par')
       split_names.append(tmp_copy)
   return split_counts, split_list, split_names



def loss(ccm, tree, rc, names, split_list, split_counts,split_names, clusters):
   m = []  #another list to hold mutated names (includes duplicates)
   for node1 in tree.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
     ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
     for node2 in tree.traverse("preorder"):
      if len(des)==0:
       break
      else:
       if node1 != node2 and node2 not in ans:
         if node2 in des:
           if node1.is_leaf():
             mutated_names1 = set([node1.name])
           else:
             mutated_names1 = set(node1.get_leaf_names(is_leaf_fn=None))
           if node2.is_leaf():
             non_mutated_names2 = set([node2.name])
           else:
             non_mutated_names2 = set(node2.get_leaf_names(is_leaf_fn=None))
           m.append( mutated_names1 - non_mutated_names2)
   mutated_names = []
   for i in m:
    if i not in mutated_names:
     mutated_names.append(i)
   for j in mutated_names:
     not_mutated_names = clusters - j
     tmp_arr = np.zeros(rc.shape[1])
     for name in names:
       if name in list(j):
         for cell in ccm[name]:
           tmp_arr[cell] = 1
     tmp_arr = tmp_arr.tolist()
     tmp_copy = tmp_arr.copy()
     if tmp_arr not in split_list:
       split_list.append(tmp_arr)
       split_counts.append(tmp_arr.count(1))
       tmp_copy.append('Loss')
       split_names.append(tmp_copy)
   return split_counts, split_list, split_names
   





def ML(ccm, tree, rc, names, one_Ls, zero_Ls):  
  # ccm : cluster cell map , 
  # rc : Read count matrix (site x cell), 
  # one_Ls : likelihood matrix of mutation occuring (sites x cells)
  # zero_Ls : likelihood of mutation not occuring (sites x cells) 
  clusters = copy.copy(names)
  clusters = set(clusters)
  split_list = []
  split_counts = []
  split_names = []

  split_counts, split_list, split_names = heterozygous(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)
  split_counts, split_list, split_names = parallel(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)
  split_counts, split_list, split_names = loss(ccm, tree, rc, names, split_list, split_counts,split_names, clusters)

  tmp_arr = np.zeros(rc.shape[1])  # The case for no mutation at all
  tmp_arr = tmp_arr.tolist()
  tmp_copy = tmp_arr.copy()
  split_list.append(tmp_arr) 
  split_counts.append(tmp_arr.count(1))
  tmp_copy.append('None')
  split_names.append(tmp_copy)

  split_list = np.array(split_list)
  split_names = np.array(split_names,dtype='object')

  indexes = np.argsort(split_counts)
  splits_sorted = split_list[indexes,:]
  split_names = split_names[indexes,:]
  new_genotypes = np.zeros((rc.shape[0], rc.shape[1]))
  L_wg = np.float128(0)
  s = time.time()
  (L_wg, new_genotypes,split_names) = chunk_ML(splits=splits_sorted, inverted_splits=1-splits_sorted, ones=np.float16(one_Ls), zeros=np.float16(zero_Ls),split_names = split_names)
  new_genotypes = new_genotypes.T
  return (new_genotypes,L_wg,split_names)
  

def ML_initialization(lc, one_Ls, zero_Ls):
# Initialisation of likelihood matrices (no mutation likelihood matrix and mutation likelihood matrix)
# Initialisation of the likelihood value
    init_L = np.float128(0)
    num = lc.shape[1]
    pos = lc.shape[0]
    lc = lc.replace(0,0.00000001)
    lc = lc.replace(1,0.99999999)
    for i in tqdm(range(pos),desc="Initializing Log likelihoods..."):
        for j in range(num):
            one_Ls[i][j] = np.float16(np.log(lc.iloc[i][j]))
            zero_Ls[i][j] = np.float16(np.log(1-lc.iloc[i][j]))
    one_Ls = one_Ls.T
    zero_Ls = zero_Ls.T
    return (init_L, one_Ls, zero_Ls)



def generate_graph_laplacian(df, nn):
# Generate Graph Laplacian from data 
    connectivity = kneighbors_graph(X=df, n_neighbors=nn, mode='connectivity')
    adjacency_matrix_s = (1/2)*(connectivity + connectivity.T)
    graph_laplacian_s = csgraph.laplacian(csgraph=adjacency_matrix_s, normed=True)#unnormalized laplacian
    graph_laplacian = graph_laplacian_s.toarray()
    return graph_laplacian



def compute_spectrum_graph_laplacian(graph_laplacian):
# Compute eigenvalues and eigenvectors and project them onto the real numbers
    eigenvals, eigenvcts = linalg.eig(graph_laplacian)
    eigenvals = np.real(eigenvals)
    eigenvcts = np.real(eigenvcts)
    return eigenvals, eigenvcts



def findBestLikelihood(iterations,resol,data_df,one_Ls,zero_Ls,best_cluster_s): 
# Find best likelihoods among a set of chosen number of clusters
    adata = sc.AnnData(data_df)
    adata
    sc.pp.neighbors(adata,use_rep = 'X')
    sc.tl.leiden(adata,resolution=resol)
    cluster_labels = np.array(adata.obs['leiden'],dtype='int')
    sil_score = silhouette_score(data_df, cluster_labels)
    ysc_s.append(cluster_labels)
    n_cluster = len(set(cluster_labels))
    print("Resolution : ",resol)
    print("Number of Clusters : ",n_cluster)
    if n_cluster in best_cluster_s:
       return (1,1,1,1,1,1)
    if n_cluster > 13:
       print("Number of clusters : ",n_cluster)
       return (0,0,0,0,0,0)
    parse_time = time.time()
    names=[]
    for i in range(n_cluster):
        name = "c"+str(i+1)
        names.append(name)
    cluster_cell_map={}
    for name in names:
        ind = int(name[1:])
        cluster_cell_map[name] = [i for i in range(len(cluster_labels)) if cluster_labels[i] == ind-1]
    initialization = None
    init_time = time.time()
    # Creating random tree topology with given cluster count
    ete_nj_init = Tree()
    ete_nj_init.populate(len(names),names)
    # Reconstruct the NJ tree given the initial results
    best_L = float("-inf")
    (new_mat, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=ete_nj_init, rc=data_df.T, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)

    if Likelihood>best_L:
        best_L = Likelihood

    n_iterations = iterations
    best_Ls = [best_L]
    stack  = [ete_nj_init]
    top_ids = set()
    top_ids.add(ete_nj_init.get_topology_id())
    best_mat = new_mat
    best_tree = ete_nj_init
    best_geno_name = split_names
    best_trees = [best_tree]
    best_geno_names = [split_names]
    for it in range(n_iterations):
        print("iteration # ",it)
        Ls = []
        ts = []
        mats = []
        genos = []
        for item_ in stack:
            tree_list = NNI.Main(in_tree=item_, N=n_cluster)
            for tree_ in tree_list:
                if tree_.get_topology_id() not in top_ids:
                    top_ids.add(tree_.get_topology_id())
                    (mat_, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=tree_, rc=lklhd_computed, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)
                    if Likelihood > best_L:
                        best_L = Likelihood
                        best_mat = mat_
                        best_tree = tree_
                        best_geno_name = split_names
                        Ls=[]
                        mats=[]
                        ts=[]
                        Ls.append(Likelihood)
                        mats.append(mat_)
                        ts.append(tree_)
                        genos.append(split_names)
                
            tree_list = SPR.Main(in_tree=item_, N=n_cluster, N_dest=n_cluster)
            for tree_ in tree_list:
                if tree_.get_topology_id() not in top_ids:
                    top_ids.add(tree_.get_topology_id())
                    (mat_, Likelihood, split_names) = ML(ccm=cluster_cell_map,tree=tree_, rc=lklhd_computed, names=names, one_Ls=one_Ls, zero_Ls=zero_Ls)
                    if Likelihood > best_L:
                        best_L = Likelihood
                        best_mat = mat_
                        best_tree = tree_
                        best_geno_name = split_names
                        Ls=[]
                        mats=[]
                        ts=[]
                        Ls.append(Likelihood)
                        mats.append(mat_)
                        ts.append(tree_)
                        genos.append(split_names)

                
        max_ = float("-inf")
        if len(Ls)!=0:
            stack = ts
        else:
            print("no more better proposed trees")
            print("terminating the search")
            break
    return best_L,best_mat,n_cluster,sil_score, best_tree, best_geno_name                       



data = args.dataset
iterations = args.treeIterations
outputPrefixName = data + "_indels_finite_sites_"
countFileName = "/home/priya/Documents/priya/AML_indels_batch_2_readcounts_gq_files/"+data+"/ReadCounts_"+data+".tsv"
qualityFileName = "/home/priya/Documents/priya/AML_indels_batch_2_readcounts_gq_files/"+data+"/GQ_"+data+".tsv"
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
min_dist_defined = args.minDist

random.seed(seed)
np.random.seed(seed)


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1)  


# # # # # # # # # # # # # # # # # # # # # # FINAL ALGORITHM # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
start_time = time.time()
df = pd.read_csv(countFileName, sep="\t",header=None)

df = df.to_numpy()
print("df looks like",df)
print("Dimensions of the input matrix: ",df.shape)


# # # # # # # # # # # # # # # # # # # # # # FILTERS # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
filters_start_time = time.time()
muts = df.shape[0]
cells = df.shape[1]-4
alt_alleles = []
retain_indices = []


for i in tqdm(range(muts),desc='Computing Alternate alleles...'):
    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
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
    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
        allele_counts = np.zeros((cells,2),dtype=int)
        ref = df[i][2]
        alt = df[i][3]
        ref_ind = 0
        counts = getAltAlleleCounts(df[i][4:])
        counts = np.array(counts.tolist())
        allele_counts = np.sum(counts,axis=0)
        ref_count = allele_counts[0]  # reference is present at index 0 in case of indels
        allele_counts[ref_ind] = -1
        if np.max(allele_counts)<=0:
            continue
        retain_indices.append(i)
        alt_alleles.append(alt)
    else:
       print("Problem occured while computing alternate alleles at pos ",i )

       


df = df[retain_indices,:]
df_copy = copy.copy(df)
alt_alleles_copy = copy.copy(alt_alleles)
print('Read count matrix shape after computing alternate alleles is: ',df.shape)



# .................................Quality Filter..............................................................................................
if isQualityFilteringOn: 
	df_gq = pd.read_csv(qualityFileName,sep='\t',header=None)
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
for i in tqdm(range(muts),desc='Computing likelihoods... '):
    alt = alt_alleles[alt_pos]
    alt_pos += 1
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)

    
    if df[i][2] in ['A','C','G','T'] and df[i][3] in ['A','C','G','T']:
       altIndex = charToIndex(alt)
       alt_depth = counts[:,altIndex]
    elif df[i][2] not in ['A','C','G','T'] or df[i][3] not in ['A','C','G','T']:
       alt_depth = counts[:,1]  # because in case of indels, the alternate depth is always present at index 1 (ref:alt)

    # benchmarking values
    alpha = (-0.000027183*coverage)+0.068567471
    beta = (0.007454388*coverage)+2.367486659
    w = (0.000548761*coverage)+0.540396786
    alpha_1 = (0.057378844*coverage)+0.669733191
    alpha_2 = (0.003233912*coverage)+0.399261625

    # tnbc values
    #alpha = (-0.000027183*coverage)+0.068567471
    #beta = (0.007454388*coverage)+2.367486659
    #w = (0.00013049539529072615*coverage) + 0.20485071465951765
    #alpha_1 = (0.0027279276611190525*coverage) + 4.563625352235754
    #alpha_2 = (0.00010204978394866682*coverage) + 0.07785441639236484
    
    l0 = betabinom.pmf(alt_depth, coverage, alpha, beta)
    l1 = (w*betabinom.pmf(alt_depth, coverage, alpha_1, alpha_1)) + ((1-w)*betabinom.pmf(alt_depth, coverage, alpha_2, alpha_2))
    l2 = betabinom.pmf(alt_depth, coverage, beta, alpha)
    
    l1 = np.where(l1<0,betabinom.pmf(alt_depth, coverage, alpha_1, alpha_1),l1)

    l0 = l0/(l0+l1+l2)
    lklhd_df[i] = np.where(coverage==0,0,1-l0)

       
print("-------After computing likelihood matrix-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Alternate Allele frequency filter.............................................................................
pos_retained = altAlleleFreqFilter(df,lklhd_df,minAltAlleleFreq)

df = df[pos_retained,:]
df_copy = df_copy[pos_retained,:]
alt_alleles = np.array(alt_alleles)
alt_alleles = alt_alleles[pos_retained]
lklhd_df = lklhd_df[pos_retained,:]
gq = gq[pos_retained,:]
print("-------After Alternate Allele Frequency filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")
print("how does df look after alt allele freq filter",df)


# .................................Variant Removal filter.............................................................................
pos_retained_after_third = variantRemovalFilter(df,varThreshold)
print('Positions retained after Variant removal filter : ',len(pos_retained_after_third))
df = df[pos_retained_after_third,:]
df_copy = df_copy[pos_retained_after_third,:]
alt_alleles = np.array(alt_alleles)
alt_alleles = alt_alleles[pos_retained_after_third]
lklhd_df = lklhd_df[pos_retained_after_third,:]
gq = gq[pos_retained_after_third,:]
print("-------After Variant removal filter-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# .................................Variant Removal mutated filter.............................................................................
final_pos_retained = []
variantRemovalMutated(df,lklhd_df,minMutFraction,alt_alleles,final_pos_retained)
end_time = time.time()
df = df[final_pos_retained,:]
df_copy = df_copy[final_pos_retained,:]
alt_alleles = alt_alleles[final_pos_retained]
gq = gq[final_pos_retained,:]
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
lklhd_df = lklhd_df[final_pos_retained,:]
pd.DataFrame(df).to_csv(outputPrefixName+'final_df.tsv',sep='\t',header=False,index=False)
print('Computing lklhd after filters....')
for i in range(df.shape[0]):
    counts = getAltAlleleCounts(df[i][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    lklhd_df[i] = np.where(coverage==0,np.nan,lklhd_df[i])

print("-------Final dimensions after all filters are applied-------")
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("\n")



# # # # # # # # # # # # # # # # # # # # # # GRAPH LAPLACIAN # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
lklhd_df = lklhd_df.T
data_df = pd.DataFrame(lklhd_df)
data_df =KNNImputer(n_neighbors=8).fit_transform(X=data_df)
pd.DataFrame(data_df).to_csv(outputPrefixName+'final_lklhds.tsv',sep='\t',header=False,index=False)

#adata = sc.AnnData(data_df)
#adata
#sc.pp.neighbors(adata,use_rep = 'X')
#sc.tl.leiden(adata,resolution = 0.75)
#cluster_labels = np.array(adata.obs['leiden'],dtype='int')
#n_clusters = len(set(cluster_labels))


start = time.time()


#print('Number of candidate clusters after applying eigen gap heuristic: ',n_clusters)
start_resolution = 0.25
end_resolution = 1.5
print('Start resolution : ',start_resolution)
print('End resolution : ',end_resolution)
lklhd_computed = pd.DataFrame(data_df.T)
n1 = lklhd_computed.shape[1]
l1 = lklhd_computed.shape[0]
zero_Ls = np.zeros((l1,n1))
one_Ls = np.zeros((l1,n1))
(uL, one_Ls, zero_Ls) = ML_initialization(lc=lklhd_computed, one_Ls=one_Ls, zero_Ls=zero_Ls)
lklhd_s = []
best_mat_s = []
best_cluster_s = []
best_tree_s = []
ysc_s = []
sil_score_s=[]
best_geno_name_s = []
inferredClusterLabelFileName = outputPrefixName+"inferred_cluster_labels"+".txt"


# # # # # # # # # # # # # # # # # # # # # # FINDING BEST LIKELIHOOD # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
for resol in [0.1,0.2,0.3,0.4,0.5,0.75,1]:
    print("At resolution :",resol)
    (best_L,best_mat, best_n,sil_score, best_tree,best_geno_name) = findBestLikelihood(iterations,resol,data_df,one_Ls,zero_Ls,best_cluster_s)
    if (best_L,best_mat, best_n,sil_score,best_tree,best_geno_name) == (1,1,1,1,1,1):
       print("This cluster is already done")
       continue
    if (best_L,best_mat, best_n,sil_score,best_tree,best_geno_name) == (0,0,0,0,0,0):
       print("Number of clusters greater! skipping!")
       continue
    print("Resolution : ",resol,"\tLikelihood :",best_L,"\tSilhouette Score :",sil_score)
    best_cluster_s.append(best_n)
    lklhd_s.append(best_L)
    best_mat_s.append(best_mat)
    sil_score_s.append(sil_score)
    best_tree_s.append(best_tree)
    best_geno_name_s.append(best_geno_name)


print('Algorithm completed in: ',time.time()-start_time,' secs')
print("Length of likelihood list : ",len(lklhd_s))
print("Length of matrices : ",len(best_mat_s))
print("Length of cluster list : ",best_cluster_s)
print("Length of silhouette scores : ",sil_score_s)
print("List of trees shape:",len(best_tree_s))
print("Best geno names shape :",len(best_geno_name_s))

print("Cluster list is :", best_cluster_s)
inferred_index = np.argmax(np.array(sil_score_s))
final_inferred_clusters = best_cluster_s[inferred_index]
best_mat = best_mat_s[inferred_index]
best_labels = ysc_s[inferred_index]
best_tree = best_tree_s[inferred_index]
best_geno = best_geno_name_s[inferred_index]

pd.DataFrame(best_geno).to_csv("Genotype configuration.tsv",sep='\t',index=False,header=False)
best_tree.render(outputPrefixName+"inferred_tree.png", w=183, units="mm")
best_tree.write(format=9, outfile=outputPrefixName+"inferred_tree.nw")

np.savetxt(inferredClusterLabelFileName,best_labels,fmt="%s", delimiter=",")
pd.DataFrame(best_mat).to_csv(inferredGenotypeFileName,index=False,header=False)
best_mat_trimmed = pd.DataFrame(best_mat).drop_duplicates()
inferred_genotypes = best_mat.T
final_indices = getRightOrderForVCF(df)
df = df[final_indices,:]
gq = gq[final_indices,:]
inferred_genotypes = inferred_genotypes[final_indices,:] 
vcfResultFileName = outputPrefixName+'outputInference.vcf'
vcfFile = open(vcfResultFileName,'w')
numCells = df.shape[1]-4
positions = df.shape[0]

print("-------Final dimensions after the complete algorithm -------")
print('Number of cells : ',df.shape[1]-4)
print('Number of genomic positions retained : ',df.shape[0])
print("Read count matrix shape : ",df.shape)
print("Genotype quality matrix shape : ",gq.shape)
print("Likelihood matrix shape : ",lklhd_df.shape)
print("Inferred genotypes matrix : ",inferred_genotypes.shape)
print("Length of alternate alleles : ",len(alt_alleles))
print("-------------------------------------------------------")
print("-------------------------------------------------------")

print('Writing inferred genotypes as VCF file...')

writeVCFFile(df,gq,vcfFile,inferred_genotypes)

print('Final inferred clusters in the tree: ',final_inferred_clusters)
print("Final inferred distinct genotypes: ",best_mat_trimmed.shape[0])

# making umap

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




with open('Final clusters inferred.txt','w') as f:
   f.write('Final inferred clusters in the tree\tFinal inferred distinct genotypes')
   f.write('\n')
   f.write(str(final_inferred_clusters) + '\t' + str(best_mat_trimmed.shape[0]))
   


