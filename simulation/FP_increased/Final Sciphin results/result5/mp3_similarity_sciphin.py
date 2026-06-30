import numpy as np
import pandas as pd
import re
from ete3 import Tree

data_no = 5
no_clusters = 10

vcf_obtained = np.loadtxt('data_fp_inc_'+str(data_no)+'_2000_'+str(no_clusters)+'.vcf',dtype='object')
numcells = vcf_obtained.shape[1]-9
genotype_arr = np.zeros((vcf_obtained.shape[0],vcf_obtained.shape[1]-9+2),dtype='object')
for i in range(vcf_obtained.shape[0]):
  genotype_arr[i][0] = vcf_obtained[i][0]
  genotype_arr[i][1] = int(vcf_obtained[i][1])
  for j in range(vcf_obtained.shape[1]-9):
    if vcf_obtained[i][9:][j].split(":")[0] == '0/1':
      genotype_arr[i][j+2] = 1
    elif vcf_obtained[i][9:][j].split(":")[0] == '0/0':
      genotype_arr[i][j+2] = 0
    else:
      print("Something wrong")

df_obtained = pd.DataFrame(genotype_arr)
df_obtained.drop([0], axis=1,inplace=True)  # dropping the 'chromosome' column
df_obtained.columns=range(numcells+1) # renaming the columns from 0 to numcells (first column being the position column)
df_obtained.sort_values(by = df_obtained.columns[0],inplace=True,ignore_index=True)  # sorting this data frame in ascending order of positions
cols = [i for i in range(1,numcells+1)]  # we need to add genotypes across all cells, hence selecting columns of cells (because first column is position column)
df_obtained['sum'] = df_obtained[cols].sum(axis=1) # finding sum across all cells
df_obtained = df_obtained[df_obtained['sum'] > 1]  # mutation should be present in at least two cells, hence choosing columns where sum is more than 1
mutated_pos = list(df_obtained[0])
all_pos_sciphin = np.array([int(i) for i in vcf_obtained[:,1]],dtype='int')
non_mutated = [i for i in all_pos_sciphin if i not in mutated_pos]
print("Number of sites mutated in lesser than two cells ", len(non_mutated))
print("Number of sites mutated in more than one cell ",len(mutated_pos))

with open('data_fp_inc_'+str(data_no)+'_2000_'+str(no_clusters)+'.gv','r') as f:
  content = f.readlines()

# The most basic filtering as per the type of gv content that sciphi holds. 
for i in range(len(content)):
  if content[i].find('label') > 0:   # where 'label' string is present i.e. node name and label information
    pattern = 'shape=box,style=filled, fillcolor=white,'  # remove the node styling part
    result_1 = re.sub(pattern,'',content[i])
    content[i] = result_1

    pattern = 'style=filled, fillcolor=grey82,'  # remove the node styling part
    result_2 = re.sub(pattern,'',content[i])
    content[i] = result_2

    pattern = r'\\n'   # remove double slash n's
    result = re.sub(pattern, ',', content[i])
    content[i] = result

    pattern = r'\n'  # remove single slash n's
    result = re.sub(pattern, '', content[i])
    content[i] = result

    pattern = r'"(\d+),' # they have given cell name before putting the labels, so that cell name needs to be removed. it looks like "cell_number,".. hence replace such an expression with empty string
    result = re.sub(pattern, '"', content[i])
    content[i] = result

    result2 = re.sub(',"]', '"]', content[i]) # since if there are some commas left at the end of a string of mutations, end that string with "]
    content[i] = result2

    result11 = re.sub('chr','',content[i]) # they gave mutation as "chr1_122" so, replace 'chr' with an empty string
    content[i] = result11

    result12 = re.sub('_',":",content[i])  # and replace '_' with colon because that is how our  gv trees are prepared
    content[i] = result12

for i in range(len(content)):
  if 'label' in content[i]:
    labels = re.findall('"(.*?)"',content[i])[0].split(',')  # extract the label string and split it into list 
    for pos in non_mutated:
      string = '1:'+str(pos) # making the mutation string of a non mutated site (means mutated in less than two cells)
      if string in labels:   
        labels.remove(string)  # remove that site which is not mutated
    index_of_colon = content[i].find('"') 
    content_after_colon = content[i][index_of_colon:] # find the content after index of colon
    string_new = ','.join(x for x in labels) # create a new label by removing the non mutated sites
    string_new = '"' + string_new + '"];\\n'
    content[i] = re.sub(content_after_colon,string_new,content[i]) # replace the old labels with new labels by replacing everything after colon with the new updated string

# This gv file now only contains those mutation labels which are mutated in more than one cell
# Now, removing the empty node

# creating a nw format tree using the node connection information from gv

with open('intermediate_nw.txt','w') as f:
  for i in range(len(content)):
    if '->' in content[i]:
      child = content[i].split('->')[0].strip()  # if any extra whitespaces are present
      parent = content[i].split('->')[1].strip() # if any extra whitespaces are present
      temp = re.sub(";",'',parent)
      parent = temp.strip()
      f.write(child)
      f.write('\t')
      f.write(parent)
      f.write('\t'+str(1.0))
      f.write('\n')

t = Tree.from_parent_child_table([line.split() for line in open("intermediate_nw.txt")] )
for node in t.traverse('preorder'):
  node.label = []  # initialise all labels with an empty list

# this node_label dictionary contains info like - node_name : [list of mutated positions in string form]
  # e.g. node_name = {1786:['1:22','1:6765']}
node_label = dict()
for i in range(len(content)):
  if "label" in content[i]:
    node_name = re.findall(r'(\d+)\[', content[i])[0]
    node_label[node_name] = re.findall('"(.*?)"',content[i])[0].split(',')

for node in t.traverse('preorder'):
  for node_key in node_label.keys():
    if node.name == node_key:
      node.label = node_label[node_key]

for node in t.traverse("preorder"):
  if node.label == ['']:
    node.delete(prevent_nondicotomic=False)

for node in t.traverse("preorder"): # arrow here refers to the arrow we use to tell the relationship between two nodes, so it basically has the
                                          # information about the children of an internal node
    node.arrow = []

# writing as gv format
for node in t.traverse("preorder"):
  if not node.is_leaf():
    children = node.get_children()
    node.arrow.extend(children)
string = ''
for node in t.traverse("preorder"):
  if not node.is_leaf():
    k=0
    for j in range(len(node.get_children())):
      string_k = str(node.name) + ' -> ' + str(node.arrow[j].name) + '\n'
      string = string + string_k
      k += 1

with open('Sciphin_fp_inc_'+str(no_clusters)+'_2000_data_'+str(data_no)+'_empty_removed.gv','w') as f:
  f.write('digraph G {\n')
  for node in t.traverse("preorder"):
    list_of_label = node.label
    delimiter = ','
    string_of_label = delimiter.join(list_of_label)
    f.write(str(node.name)+' [label="'+string_of_label+'"];\n')
  f.write(string)
  f.write("}")
print(" nw tree written as gv format")
