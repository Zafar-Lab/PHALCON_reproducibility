import numpy as np
import pandas as pd
import pickle
import random

data = '2'
seed = 151161
sample_name = "Compass_input_data_2000_cpn_025_"+ str(data)
region = dict()
genome_length = 0


read_count_matrix = pd.read_csv('sc_2000_readCounts_'+ data +'.tsv',sep="\t",header=None)
read_count_np = np.array(read_count_matrix, dtype= 'object')

cells = read_count_np.shape[1] - 4
sites = read_count_np.shape[0]

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


number_of_genes = random.randint(14,18)
genomic_fragments = list(split(range(sites), number_of_genes))
print(genomic_fragments)


with open('genomic_fragments.pkl','wb') as f:
    pickle.dump(genomic_fragments, f)



def getAltAlleleCounts(arr):
    return np.array(list(map(int,arr.strip().split(","))))


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1) 

for site in range(sites):
    counts = getAltAlleleCounts(read_count_np[site][4:])
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    read_count_np[site][4:] = coverage


print(read_count_np)
read_count = pd.DataFrame(read_count_np)
read_count.drop([0,2,3],axis=1, inplace=True)
read_count.to_csv("Total_depth_read_count.tsv", sep = '\t', header=False, index=False)


read_count = pd.read_csv("Total_depth_read_count.tsv", sep = '\t', header=None)
gene_split = list()
for gene in genomic_fragments:
    print("gene",gene)
    gene_split.append(list(split(gene,25)))

with open('gene_split.pkl', 'wb') as f:
    pickle.dump(gene_split, f)

read_count = pd.read_csv("Total_depth_read_count.tsv", sep = '\t', header=None)
with open('gene_split.pkl','rb') as f:
    gene_split = pickle.load(f)

with open('genomic_fragments.pkl','rb') as f:
    genomic_fragment = pickle.load(f)



average = []
final_df = []
region_sum = []
for i in range(len(genomic_fragment)):
    one_row = []
    average = []
    for gene_frag in gene_split[i]:
        slice_df = read_count.iloc[list(gene_frag)]
        average.append(slice_df.mean())
    averages_df = pd.DataFrame(average) 
    region = str(1) + "_Region" + str(i)
    one_row.append(region)
    one_row.extend(list(map(int, list(averages_df.sum(numeric_only=True, axis=0)[1:2001]))))
    final_df.append(one_row)

final_df = pd.DataFrame(final_df)
print(final_df.isnull().values.any())

print("Final regions dataframe:", final_df)

final_df.to_csv(sample_name+"_regions.csv",header=False, index=False)

