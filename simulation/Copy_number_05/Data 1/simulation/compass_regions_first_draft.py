from ete3 import Tree
from scipy.stats import beta
from scipy.stats import multinomial
import random
import numpy as np
import pandas as pd
import copy
from sklearn.utils import shuffle
import pickle

data = '1'
sample_name = "Compass_input_data_2000_cpn_05_"+ str(data)
# selecting regions from the genome of length 40000
region = dict()
genome_length = 0


read_count_matrix = pd.read_csv('sc_2000_readCounts_'+ data +'.tsv',sep="\t",header=None)
read_count_np = np.array(read_count_matrix, dtype= 'object')

cells = read_count_np.shape[1] - 4
sites = read_count_np.shape[0]

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


genomic_fragments = list(split(range(sites), 15))
print(genomic_fragments)

# finding the total read depth for every site and every

def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
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
        #print("gene frag:",gene_frag)
        slice_df = read_count.iloc[list(gene_frag)]
        #print("slice df",slice_df)
        average.append(slice_df.mean())
    averages_df = pd.DataFrame(average)
    #print("averages df",averages_df)
    
    region = str(i) + "_REGION" + str(i)
    one_row.append(region)

    one_row.extend(list(map(int, list(averages_df.sum(numeric_only=True, axis=0)[1:2001]))))

    final_df.append(one_row)

final_df = pd.DataFrame(final_df)

print("Final regions dataframe:", final_df)
#averages_df = pd.DataFrame(average)
final_df.to_csv(sample_name+"_regions.csv",sep=',',header=False, index=False)



#average = []
#for small_region in gene_split:
 #   slice_df = read_count.iloc[list(small_region)]
  #  average.append(slice_df.mean())

#averages_df = pd.DataFrame(average)
#print(averages_df)
#    
#for gene_region in gene_split:
 #       ind_list = gene_region
  #      print("ind list",ind_list)
   #     read_count_region = read_count.iloc(list(ind_list))
     #   print("read count region",read_count_region)
      #  print("mean",read_count_region.mean(read_count_region, axis = 0))







'''
for i in range(20):
    random_length = random.randint(2500, 3000)
    if genome_length > 40000 - 2500:
        random_length = 40000 - genome_length
    region[i] = list(genome_length + np.array(range(0,random_length+1)))
    genome_length += random_length

    if genome_length == 40000:
        break


# finding the total read depth for every site and every

def getAltAlleleCounts(arr):
# Converts read count strings (e.g. "0,1,2,1") to numpy array (e.g. [0,1,2,1])
    return np.array(list(map(int,arr.strip().split(","))))


getAltAlleleCounts =np.frompyfunc(getAltAlleleCounts, 1, 1) 

for site in range(sites):
    
    counts = getAltAlleleCounts(read_count_np[site][4:])
    
    counts = np.array(counts.tolist())
    coverage = np.sum(counts,axis=1)
    read_count_np[site][4:] = coverage


    '''


''' 
for cell in range(cells):
    read_depth = 0
    for site in range(sites):
        counts = read_count_np[site][cell+4]        
        counts =counts.split(",")
        counts = list(counts)
        counts = [int(counts[i]) for i in range(len(counts))]
    
        read_depth = sum(counts)
        read_count_np[site][cell+4] = read_depth

print(read_count_np)

for region_key in region.keys():
    print("region keys:",region_key)
    for site in region[region_key]:
        print(site)
        site_list = list(split(range(site), 3))
        print(site_list)

        #if site in region[region_key]:
         #   read_count_np[site][0] = region_key
          #  read_count_np[site][1] = str(region_key) + "_Region" + str(region_key)



read_count = pd.DataFrame(read_count_np)

read_count.drop([2,3],axis=1, inplace=True)


print(read_count)


read_count = read_count.groupby(1).max()
print(read_count.columns)
read_count_region = read_count.drop([0],axis=1)



print(read_count_region)
read_count_region.to_csv("Region_file_temp.csv", sep = ',', header=False)
'''

'''
region_file = []
for region_no, region_sites in region.items():
    region_row = []
    string = str(region_no) + "_REGION" + str(region_no)
    region_row.append(string)

    for cell in range(cells):
        read_depth = 0
        for site in range(sites):
            if site >= region_sites[0] and site <= region_sites[1]:
                read_depth += read_count_np[site][cell+4]
            region_row.append(read_depth)
    region_file.append(region_row)

region_file = pd.DataFrame(region_file)
region_file.to_csv("Region_file_temp.csv", sep = ',', header=False, index=False)
                            
'''


