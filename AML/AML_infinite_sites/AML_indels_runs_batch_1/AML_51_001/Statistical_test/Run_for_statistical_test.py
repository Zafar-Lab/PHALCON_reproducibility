import os
import ast
import pandas as pd
import argparse

argParser = argparse.ArgumentParser(prog='PROG')

argParser.add_argument('-D', '--Dataset', type=str)

args = argParser.parse_args()

data = args.Dataset
df_for_sites = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/AML_indels_nonbulk_evident_sites_dbsnp_nonzero_statistical_BF_test.tsv",sep='\t')
df_for_sites_of_a_sample = df_for_sites[df_for_sites['Sample_id'] == data]
index_for_the_site = df_for_sites_of_a_sample.index[0]
print(df_for_sites_of_a_sample.index)
print("this\n",df_for_sites_of_a_sample['Sites that are not present in bulk'])
list_bulk_evidence_not_found = ast.literal_eval(df_for_sites_of_a_sample['Sites that are not present in bulk'][index_for_the_site])

print("List of sites which are not found in bulk evidence:",list_bulk_evidence_not_found)

with open('Log BF test.txt','w') as f:
   f.close()

with open('Log BF test.txt','a+') as f:
   f.write("Data\tSite number\tPHALCON_inferred_log_lklhd\tFP_inferred_log_lklhd\n")

for site in list_bulk_evidence_not_found:
  print(site)
  os.system("python Algorithm_phalcon.py -S {} -d {}".format(site,data))
  os.system("python Algorithm_sciphin.py -S {} -d {}".format(site,data))

