import re
import os
import pandas as pd

'''
This python file will go to each AML directory and run the files that are stated
'''

folder = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

for aml_file_name in sub_folders:
    print("File working on: ",aml_file_name)
    os.chdir(folder+aml_file_name+"/")
    os.system("python3 GV_tree_cluster_mismatch_dbsnp_gene_tree_exonic.py")
