from ete3 import Tree
from scipy.stats import beta
from scipy.stats import multinomial
import random
import numpy as np
import pandas as pd
import copy
from sklearn.utils import shuffle

seed = 152322
numCells = 2000
numClones = 10
alpha_cp = 1
beta_cp = 1
total_pos = 40000     # Total number of positions
mutated_pos = 0.00125    #Fraction of mutated positions out of total positions
m_pos = total_pos * mutated_pos  # number of mutated positions
print("Number of mutated positions :",m_pos)
phi = ['None'] * numClones
for i in range(numClones):
  r = beta.rvs(alpha_cp, beta_cp, size = 1)[0]
  phi[i] = r
phi = [phi[i]/sum(phi) for i in range(len(phi))]  # normalising the beta prior to sum up to one

cell_in_clone = list(np.random.multinomial(numCells, phi, size=1)[0])  # sampling number of cells from each clone based on their probability vector

# assigning cells to clones based on the the above list
cluster_cell_map={}
cells = list(range(1, numCells+1))
count = 0
names = []
for j in cell_in_clone:
  name = "C" + str(count+1)
  cluster_cell_map [name] = []
  random_sample = random.sample(cells, j)
  for i in random_sample:
    cells.remove(i)
  cluster_cell_map[name] = random_sample
  count += 1
  names.append(name)
t = Tree()  # intialize a tree with 10 clones (numClones number of clones)
t.populate(numClones,names)  # this is our TRUE TREE
t.render("true tree.png",w=183,units='mm' )
t.write(format=9, outfile="data_ado_inc_"+str(numCells)+"_"+str(numClones)+"_inferred_tree.nw")

print("Length of ccm : (should be equal to the number of clusters) :",len(cluster_cell_map) )

# all cases of parallel mutation
m = []  #list to hold mutated names (includes duplicates)
for node1 in t.traverse("preorder"):
    des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
    ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
    for node2 in t.traverse("preorder"):
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
mutated_names_par = []
for i in m:
    if i not in mutated_names_par:
      mutated_names_par.append(i)

# all cases of back mutation
m = []
for node1 in t.traverse("preorder"):
  des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
  ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
  for node2 in t.traverse("preorder"):
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
mutated_names_back = []
for i in m:
  if i not in mutated_names_back:
    mutated_names_back.append(i)

# all cases of heterozygous mutation
m = []
for node in t.traverse("preorder"):
  if node.is_leaf():
    mutated_names = set([node.name])
  else:
    mutated_names = set(node.get_leaf_names(is_leaf_fn=None))
  m.append(mutated_names)
mutated_names_het = []
for i in m:
  if i not in mutated_names_het:
    mutated_names_het.append(i)

    # Adding code for back mutations of clonal cases only
m = []
first_node = 0
for node1 in t.traverse("preorder"):
  first_node += 1
  des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
  ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
  for node2 in t.traverse("preorder"):
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
  if first_node == 1:
    break
mutated_names_back_clonal = []
for i in m:
  if i not in mutated_names_back_clonal:
    mutated_names_back_clonal.append(i)

def split_list_fn(mutated_names):
  split_list = []
  for j in mutated_names:
    tmp_arr = np.zeros(numCells)
    for name in names:
      if name in list(j):
        for cell in cluster_cell_map[name]:
          tmp_arr[cell-1] = 1
    tmp_arr = tmp_arr.tolist()
    split_list.append(tmp_arr)
  return split_list

split_list_par = split_list_fn(mutated_names_par)
split_list_back = split_list_fn(mutated_names_back)
split_list_het = split_list_fn(mutated_names_het)
split_list_clonal_back = split_list_fn(mutated_names_back_clonal)

split_list = []
for i in split_list_back:
  if i not in split_list:
    split_list.append(i)
for i in split_list_het:
  if i not in split_list:
    split_list.append(i)
for i in split_list_par:
  if i not in split_list:
    split_list.append(i)
for i in split_list_clonal_back:
  if i not in split_list:
    split_list.append(i)

random.seed(seed)
TG_mut = [None] * int(m_pos) # true genotype of mutated positions
print("Length of TG_mut (should be equal to the true number of mutations)",len(TG_mut))
geno_dict = {}
TG_type = ['None'] * int(m_pos)

for i in range(int(m_pos)):   # at remaining positions putting back, parallel and heterozygous bases on some probability
  if random.random() < 0.3:  # 30% clonal : out of that, 95% clonal and 5% back
    if random.random() < 0.05:
      TG_mut[i] = random.sample(split_list_clonal_back,1)[0]  # Putting back mutations that can happen in clonal mutation case
      TG_mut_list = copy.copy(TG_mut[i])
      TG_mut_list.append('Back')
      geno_dict[i] = TG_mut_list
    else:
      TG_mut[i] = list(np.ones((1,numCells))[0]) # Putting clonal mutations at 95% of the sites
      TG_mut_list = copy.copy(TG_mut[i])
      TG_mut_list.append('Het')
      geno_dict[i] = TG_mut_list
  else:
    if random.random() < 0.1:
      if random.random() < 0.5:
        TG_mut[i] = random.sample(split_list_back,1)[0]
        TG_mut_list = copy.copy(TG_mut[i])
        TG_mut_list.append('Back')
        geno_dict[i] = TG_mut_list
      else:
        TG_mut[i] = random.sample(split_list_par,1)[0]
        TG_mut_list = copy.copy(TG_mut[i])
        TG_mut_list.append('Par')
        geno_dict[i] = TG_mut_list
    else:
      TG_mut[i] = random.sample(split_list_het,1)[0]
      TG_mut_list = copy.copy(TG_mut[i])
      TG_mut_list.append('Het')
      geno_dict[i] = TG_mut_list
print("Length of geno_dict (should be equal to number of mutated positions)",len(geno_dict))
# checking the number of clonal mutations
c=0
for i in TG_mut:
  if i == list(np.ones((1,numCells))[0]):
    c+=1
print("Number of clonal mutations :",c)
genotyped_df = pd.DataFrame(geno_dict).T
genotyped_df = shuffle(genotyped_df)
genotyped_df = genotyped_df.reset_index(drop=True)
genotyped_df.to_csv("Genotyped_dataframe_"+ str(numCells) + "_" + str(int(m_pos)) + ".csv",header=None, index=None)
geno_df = pd.read_csv("Genotyped_dataframe_"+ str(numCells) + "_" + str(int(m_pos)) + ".csv",header=None)
TG_vectors = geno_df.drop(geno_df.columns[-1], axis=1)
TG_vectors = np.array(TG_vectors,dtype='int')
pd.DataFrame(TG_vectors).to_csv("Genotype_data_"+ str(numCells) + "_" + str(int(m_pos)) + ".tsv", sep = '\t',header = None)
print("Number of mutated positions: ",TG_vectors.shape[0])
print("TG_vectors shape (should be equal to the number of cells): ",TG_vectors.shape[1])   # should be equal to number of cells
ccm_list = list(cluster_cell_map.items())
ccm_dict = dict()
n = 0
for i in ccm_list:
  j = list(i)
  ccm_dict[n] = j[1]
  n = n+1
clones = []
for i in range(1,numCells+1):
  v = [t for t in ccm_dict if i in ccm_dict[t]]
  clones.append(v[0])
print("Length of clones (should be equal to number of cells) :",len(clones))
path_for_ccm = "ccm_" + str(numCells) + "_" + str(int(m_pos)) + ".txt"
with open(path_for_ccm, 'w') as f:
  for i in range(1,numCells+1):
    f.write(str(i-1)+' ')
  f.write('\n')
  for i in range(1,numCells+1):
    f.write(str(clones[i-1])+' ')