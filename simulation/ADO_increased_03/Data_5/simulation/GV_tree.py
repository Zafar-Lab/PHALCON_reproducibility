import pandas as pd
import numpy as np
from ete3 import Tree
from graphviz import Source
import re

# CHANGE THE TREE FROM HERE
tree_nw = Tree("inferred_format.nw")
data = 5

print(tree_nw)



def findPar(tree_nw, split_list):
   for node1 in tree_nw.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  # descendants of current node
     ans = list( h for h in node1.get_ancestors() )  # anscestors of current node
     for node2 in tree_nw.traverse("preorder"):
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
           m = mutated_names1.union( mutated_names2)
           if set((np.array(list(m),dtype='int'))) == set([i for i,val in enumerate(split_list) if val==1]):
            return node1, node2
def findLoss(tree_nw, split_list):
   for node1 in tree_nw.traverse("preorder"):
     des = list( j for j in node1.iter_descendants() )  #list of descendants of current node
     ans = list( h for h in node1.get_ancestors() )  #list of anscestors of current node
     for node2 in tree_nw.traverse("preorder"):
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
           m = mutated_names1 - non_mutated_names2
           if set((np.array(list(m),dtype='int'))) == set([i for i,val in enumerate(split_list) if val==1]):
            return node1, node2
def findNode(tree_nw, split_list, config):
  if config == 'Het' or config == 'Clonal':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    if len(mutated) == 1:
      return tree_nw & str(mutated[0])
    return tree_nw.get_common_ancestor([str(i) for i in mutated])
  elif config == 'Par':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    node1, node2 = findPar(tree_nw, split_list)
    return node1, node2
  elif config == 'Loss' or config == 'Back':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    node1, node2 = findLoss(tree_nw, split_list)
    return node1, node2
  else:
    print("Problem")
    return 0
# for inferred files, we do this because the cluster information here is stored in a format where there is no line for cell number, whereas
# in case of original file, the first line tells the cell number and seconf line contains the custer information and hence in that case
# we need to go to a second level indexing
config = pd.read_csv('Genotyped_dataframe_2000_50.csv',sep=',',header=None) # Genotype configuration file
orig_pos = pd.read_csv('sc_2000_inclusion_output.vcf',sep='\t',header=None) # final data frame for finding the sites which are mutated
config.insert(loc = 0,column = 'site',value = orig_pos[1]) # putting the site information at the first column
config2 = config.T.drop_duplicates().T # selecting the first occurence of a different genotype across all cells, i.e. finding unique genotypes  which will readily convert the cellular dataframe into clonal one
cluster_labels = np.loadtxt('ccm_2000_50.txt',dtype='int') # now finding the labels corresponding to each unique cell
cluster_labels = np.delete(cluster_labels, (0), axis=0)
cluster_labels = cluster_labels[0]
columns = dict()
for i in config2.columns[1:-1]:
  columns[i] = cluster_labels[i]
config2.rename(columns = columns, inplace = True)  # replace the cell names with the clonal information obtained just above

'''Need to change it according to the number of clusters you have'''
# CHANGE THE CODE BELOW AS PER THE NUMBER OF CLUSTERS
cols = []
cols.append('site')
for i in range(len(set(cluster_labels))):
  cols.append(i)
cols.append(len(cluster_labels))  # reorder the columns so that there is no confusion 
config2 = config2[cols]
#config2 : clonal level configuration ready
clonal_info = np.array(config2,dtype='object')  # into numpy format for easy element access

i=0
for node in tree_nw.traverse("preorder"): # naming the nodes because by default they are just empty strings
    node.temp = i+len(set(cluster_labels))
    i+=1
for node in tree_nw.traverse("preorder"): # arrow here refers to the arrow we use to tell the relationship between two nodes, so it basically has the 
                                          # information about the children of an internal node
    node.arrow = []
for node in tree_nw.traverse("preorder"): # initialise the labels with an empty list
    node.label = []
for i in range(clonal_info.shape[0]): # every site individually
    site = clonal_info[i][0]
    config = clonal_info[i][-1]
    split_list = clonal_info[i][1:-1]
    if config == 'Clonal' or config == 'Het':
      node = findNode(tree_nw, split_list, config)
      node.label.append('1:'+str(site))
    elif config == 'Par' or config == 'Back' or config == 'Loss':
      node1,node2 = findNode(tree_nw, split_list, config)
      node1.label.append('1:'+str(site))
      if config == 'Back' or config == 'Loss':
        node2.label.append('1:'+str(site)+'-')
      else:
        node2.label.append('1:'+str(site))
    else:
      print("Fault")
for node in tree_nw.traverse("preorder"):
  print(node.label)


for node in tree_nw.traverse("preorder"):
  if not node.label:
    node.delete(prevent_nondicotomic=False)


for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    node.arrow.extend(children)
string = ''
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    k=0
    for j in range(len(node.get_children())):
      string_k = str(node.temp-len(set(cluster_labels))) + ' -> ' + str(node.arrow[j].temp-len(set(cluster_labels))) + '\n'
      string = string + string_k
      k += 1



with open('True_tree_ado_inc_2000_data_'+str(data)+'.gv','w') as f:
  f.write('digraph G {\n')
  for node in tree_nw.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.temp-len(set(cluster_labels)))+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")
