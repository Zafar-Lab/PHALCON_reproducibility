import pandas as pd
import numpy as np
from ete3 import Tree
from graphviz import Source
import re
import os


data=os.path.basename(os.getcwd())

with open(data+'_indels_inferred_tree.nw','r') as f:
        string = f.read()


print(string)

print(string)
string = re.sub('c','',string)
print(string)
with open('inferred_format.nw','w') as f:
        f.write(string)


def decrease_numbers(text):
    def replace(match):
        num = match.group(0)
        return str(int(num)-1)

    pattern = r'\b\d+\b'
    return re.sub(pattern, replace, text)

# Read the text file
file_path = 'inferred_format.nw'  # Update this with the path to your text file
with open(file_path, 'r') as file:
    original_text = file.read()

# Decrease numbers in the text
modified_text = decrease_numbers(original_text)

# Write the modified text back to the file
with open(file_path, 'w') as file:
    file.write(modified_text)
    file.write(";")

with open('inferred_format.nw','r') as f:
        content = f.read()
        print(content)


        



# CHANGE THE TREE FROM HERE
# Use cluster information from outside the folder
tree_nw = Tree("inferred_format.nw")




def findNode(tree_nw, split_list, config):
  if config == 'Het' or config == 'Clonal':
    mutated = [i for i,val in enumerate(split_list) if val==1]
    if len(mutated) == 1:
      return tree_nw & str(mutated[0])
    return tree_nw.get_common_ancestor([str(i) for i in mutated])
  #elif config == 'Par':
   # mutated = [i for i,val in enumerate(split_list) if val==1]
    #node1, node2 = findPar(tree_nw, split_list)
    #return node1, node2
  #elif config == 'Loss' or config == 'Back':
   # mutated = [i for i,val in enumerate(split_list) if val==1]
    #node1, node2 = findLoss(tree_nw, split_list)
    #return node1, node2
  else:
    print("Problem")
    return 0
# for inferred files, we do this because the cluster information here is stored in a format where there is no line for cell number, whereas
# in case of original file, the first line tells the cell number and seconf line contains the custer information and hence in that case
# we need to go to a second level indexing

dbsnp_non_zero_avinput_file_name = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_2/'+data+'/'+data+'_indels_dbsnp_non_zero.avinput.variant_function'
final_variants_df = pd.read_csv(dbsnp_non_zero_avinput_file_name,sep='\t',header=None)
# final_variants = final_variants[3] -> this is okay in case of insertions and SNVs but not in deletions


# this extra piece of code is added in order to include deletions also, because annovar handles deletions a little differently
# annovar is like if at position 123 -> CT changes to T it means that at position 124 T changes to - i.e. a deletion has happened of T
# so, we check in annovar output file, if annovar_file column 6 (0-based index), we have a - , we convert that variant position into position-1
# because our read count file is made in this way ( position 123 -> CT changes to T ) rather than (position 124 T changes to -)
final_variants = []
final_variants_df = final_variants_df.to_numpy()
for i in final_variants_df:
   if i[6] == '-':
      final_variants.append(i[3]-1)
   else:
      final_variants.append(i[3])
   
   

orig_pos = pd.read_csv(data+'_indels_final_df.tsv',sep='\t',header=None) # this is the dataframe which is obtained after post processing and everything
orig_pos = orig_pos[orig_pos[1].isin(final_variants)]

config = pd.read_csv('Genotype configuration.tsv',sep='\t',header=None) # Genotype configuration file
before_post_process_file = pd.read_csv(data+'_indels_final_df.tsv',sep='\t',header=None) # final data frame for finding the sites which are mutated
# The above dataframe is the one which is obtained just after phalcon is run on the dataset
numcells = config.shape[1]-1
print("Number of cells :",numcells)

config.insert(loc = 0,column = 'chr', value = before_post_process_file[0])  # putting chromosome information at first column
config.insert(loc = 1,column = 'site',value = before_post_process_file[1]) # putting the site information at the second column

config = config[config['site'].isin(orig_pos[1])]

print("Shape of configuration file:",config.shape)
cols = [i for i in range(0,numcells)]

config['sum'] = config[cols].sum(axis=1) # sum across all cells
config = config[config['sum'] != 0]  # if sum is more than 0, keep those sites i.e. if more than 0 cell is mutated i.e. at least one cell should be mutated
config.reset_index(inplace=True,drop=True)
config  = config.drop(['sum'],axis=1)

config2 = config.T.drop_duplicates().T # selecting the first occurence of a different genotype across all cells, i.e. finding unique genotypes  which will readily convert the cellular dataframe into clonal one

cluster_labels = list(np.loadtxt(data+'_indels_inferred_cluster_labels.txt',dtype='int')) # now finding the labels corresponding to each unique cell
columns = dict()
for i in config2.columns[2:-1]:
  columns[i] = cluster_labels[i]

no_of_clusters = len(set(cluster_labels))
print("Number of clusters:",no_of_clusters)

missing = 0

for cluster_no in range(no_of_clusters):
  if cluster_no not in columns.values():
    cell_instance = cluster_labels.index(cluster_no) # index at which the first instance of that cluster appears in the crowd of all cells
    config2.insert(loc = config2.shape[1]-1,column = cell_instance,value = config[cell_instance])
    columns[cell_instance] = cluster_no
    missing += 1
    # find an instance of that cluster number in the cell configuration and add in the config
print("No of Missing columns : ",missing)


config2.rename(columns = columns, inplace = True)  # replace the cell names with the clonal information obtained just above



cols = []
cols.append('chr')
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
    chr = clonal_info[i][0][3:]
    site = clonal_info[i][1]
    config = clonal_info[i][-1]
    split_list = clonal_info[i][2:-1]
    if config == 'Clonal' or config == 'Het':
      node = findNode(tree_nw, split_list, config)
      node.label.append(str(chr) + ":" + str(site))
    #elif config == 'Par' or config == 'Back' or config == 'Loss':
     # node1,node2 = findNode(tree_nw, split_list, config)
      #node1.label.append(str(chr) + ":" + str(site))
      #if config == 'Back' or config == 'Loss':
       # node2.label.append(str(chr) + ":" + str(site) + '-')
      #else:
      #  node2.label.append(str(chr) + ":" + str(site))
    else:
      print("Fault")
for node in tree_nw.traverse("preorder"):
  print(node.label)

#########################################################This is extra code to not remove the empty nodes in phalcon inferred tree
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

with open(data+'_indels_dbsnp_nz_inferred_tree.gv','w') as f:
  f.write('digraph G {\n')
  for node in tree_nw.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.temp-len(set(cluster_labels)))+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")
########################################################################################################################################33

for node in tree_nw.traverse("preorder"):
  if not node.label:
    node.delete(prevent_nondicotomic=False)



for node in tree_nw.traverse("preorder"): # arrow here refers to the arrow we use to tell the relationship between two nodes, so it basically has the 
                                          # information about the children of an internal node
    node.arrow = []


for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    node.arrow.extend(children)

evolution_pattern='Linear'
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    if len(children) > 1:
       evolution_pattern = 'Branching'
       break
    
    


string = ''
for node in tree_nw.traverse("preorder"):
  if not node.is_leaf():
    k=0
    for j in range(len(node.get_children())):
      string_k = str(node.temp) + ' -> ' + str(node.arrow[j].temp) + '\n'
      string = string + string_k
      k += 1


all_pos = []
for node in tree_nw.traverse("preorder"):
  all_pos.extend(node.label)

print("Number of labels ",len(all_pos))

with open(data+'_indels_dbsnp_nz_nonempty_inferred_tree.gv','w') as f:
  f.write('digraph G {\n')
  for node in tree_nw.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.temp)+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")

with open(data+'_indels_evolution_pattern.txt','w') as f:
  f.write(data+'\t'+evolution_pattern)

def gene_tree_non_empty(data): # this function will convert all mutations named as "1:29997234" into their gene names
  gene_name = pd.read_csv(data+'_indels_dbsnp_non_zero.avinput.variant_function',sep='\t',header=None)
  with open(data+'_indels_dbsnp_nz_nonempty_inferred_tree.gv','r') as f:
     content = f.readlines()
  gene_name = np.array(gene_name,dtype='object')
  gene_pos = dict()
  for i in range(gene_name.shape[0]):
    chr = gene_name[i][2][3:]
    if gene_name[i][6] == '-':
       site = gene_name[i][3]-1
    else:
       site = gene_name[i][3]
    key = str(chr) + ":" + str(site)  # creating the mutation string eg: "1:204710"
    value = gene_name[i][1]
    gene_pos[key] = value
  for i in range(len(content)):
    if 'label' in content[i]:
      for position in gene_pos.keys():
        gene_name = gene_pos[position]
        pattern = r":(.+)" # extract everything after the colon, this is regex expression for that
        match = re.search(pattern, position)
        site_referred_to = match.group(1)
        tmp = re.sub(str(position),str(gene_name)+"_"+str(site_referred_to),content[i])
        content[i] = tmp
  for i in range(len(content)):
    print(content[i])

  with open(data+'_indels_dbsnp_nz_nonempty_gene_tree.gv','w') as f:
    for i in range(len(content)):
      f.write(content[i])
  print("nw tree written as gv format")

  s = Source.from_file(data+'_indels_dbsnp_nz_nonempty_gene_tree.gv')
  s.render()


def gene_tree(data): # this function will convert all mutations named as "1:29997234" into their gene names
  gene_name = pd.read_csv(data+'_indels_dbsnp_non_zero.avinput.variant_function',sep='\t',header=None)
  with open(data+'_indels_dbsnp_nz_inferred_tree.gv','r') as f:
     content = f.readlines()
  gene_name = np.array(gene_name,dtype='object')
  gene_pos = dict()
  for i in range(gene_name.shape[0]):
    chr = gene_name[i][2][3:]
    if gene_name[i][6] == '-':
       site = gene_name[i][3]-1
    else:
       site = gene_name[i][3]
    key = str(chr) + ":" + str(site)  # creating the mutation string eg: "1:204710"
    value = gene_name[i][1]
    gene_pos[key] = value

  for i in range(len(content)):
    if 'label' in content[i]:
      for position in gene_pos.keys():
        gene_name = gene_pos[position]
        pattern = r":(.+)" # extract everything after the colon, this is regex expression for that
        match = re.search(pattern, position)
        site_referred_to = match.group(1)
        tmp = re.sub(str(position),str(gene_name)+"_"+str(site_referred_to),content[i])
        content[i] = tmp
  for i in range(len(content)):
    print(content[i])

  with open(data+'_indels_dbsnp_nz_gene_tree.gv','w') as f:
    for i in range(len(content)):
      f.write(content[i])
  print("nw tree written as gv format")

  s = Source.from_file(data+'_indels_dbsnp_nz_gene_tree.gv')
  s.render()

gene_tree_non_empty(data)
gene_tree(data)