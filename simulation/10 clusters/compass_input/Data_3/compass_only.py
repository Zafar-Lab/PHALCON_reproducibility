import pandas as pd
import numpy as np
import re

data_no = '3'
no_clusters = '10'
orig_df = pd.read_csv('/home/priya/Downloads/Final Stage/Benchmarking/'+no_clusters+' clusters/Data '+data_no+'/simulation/sc_2000_readCounts_'+data_no+'.tsv',header=None, sep='\t')
gv_modified_compass = '/home/priya/Downloads/Final Stage/COMPASS/2000_cells_'+no_clusters+'_clusters/Data_'+data_no+'/Set 1/Compass_2000_'+no_clusters+'_'+data_no+'_tree_modified.gv'  
# this is the gv file which was inferred by compass (all beautifications from this file have been removed)


with open(gv_modified_compass) as f:
  content = f.readlines()

# since we do not need node attachment information, we only keep the textual content 
  # which has the label information
new_content = []
for i in content:
  if 'label' in i:
    new_content.append(i)


# remove instances such as "0[label = <"
# remove ">];"
# remove "1:"
# remove "\\n"
# after this, we will be left with only site numbers in text form in a list
for i in range(len(new_content)):
  line1 = re.sub('(\d+)\[label=<','',new_content[i])
  new_content[i] = line1
  line2 = re.sub('>];','',new_content[i])
  new_content[i] = line2
  line3 = re.sub("1:",'',new_content[i])
  new_content[i] = line3
  line4 = re.sub("\\n",'',new_content[i])
  new_content[i] = line4

print(new_content)

compass_sites = []
for i in new_content:
  if i==' ' or '':  # if there is no label, some instances empty string might also come into picture
    continue
  else:
    sites = i.split(',')
    sites = [int(site) for site in sites]
    compass_sites.extend(sites)


print("True number of sites : ",len(compass_sites))
df = orig_df[orig_df[1].isin(compass_sites)]
df.to_csv("sc_2000_readCounts_"+no_clusters+"_"+data_no+"_compass_only.tsv",sep='\t',header=False,index=False)
# we need to convert it into a specific format as is given in their git

snv = pd.read_csv("sc_2000_readCounts_"+no_clusters+"_"+data_no+"_compass_only.tsv",sep='\t',header=None)
snv_np = np.array(snv,dtype='object')
base_to_num = {"A":0,"C":1,"G":2,"T":3}
data = []
for i in range(snv_np.shape[0]):
  chr = snv_np[i][0][3]
  ref = snv_np[i][2]
  alt = snv_np[i][3]
  name = str(chr)+":"+str(snv_np[i][1])
  freq = 0
  region = "Region1"
  row = [chr, ref, alt, region, name, freq]
  for j in range(snv_np.shape[1]-4):
    count = list(map(int,snv_np[i][4+j].strip().split(",")))
    ref_count = count[base_to_num[ref]]
    alt_count = count[base_to_num[alt]]
    row.append(str(ref_count)+":"+str(alt_count))
  data.append(row)
compass_input = pd.DataFrame(data)
columns = {0:"CHR",1:"REF",2:"ALT",3:"REGION",4:"NAME",5:"FREQ"}
# the cells should be named from 0 to numcells-1 and hence the code below:
keys = [i for i in range(6,compass_input.shape[1])]
values = [i for i in range(snv_np.shape[1]-4)]
for i in range(len(keys)):
    columns[keys[i]] = values[i]

# make the input exactly the same as compass has asked for
compass_input.rename(columns = columns,inplace=True)

# save this new format 
compass_input.to_csv("Compass_only_input_data_2000_"+no_clusters+"_"+data_no+"_variants.csv",sep=',',index=False,header=True)

# maybe see what you are sending as input
input = pd.read_csv("Compass_only_input_data_2000_"+no_clusters+"_"+data_no+"_variants.csv",sep=",")
print(input)