import re
import os
import pandas as pd
import random
'''
This python file will go to each AML directory and run the files that are stated
'''

folder = '/home/priya/Documents/priya/AML_indels_runs_batch_2/'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
sub_folders=['AML_94_001', 'AML_77_001', 'AML_67_001_PDX1', 'AML_98_001', 'AML_10_001', 'AML_99_002', 'AML_80_001', 'AML_24_001', 'AML_84_001', 'AML_108_001', 'AML_109_001', 'AML_106_001', 'AML_88_001_07', 'AML_100_001']
for aml_file_name in sub_folders:
    random_number=random.randint(1,1000)
    print("File working on: ",aml_file_name)
    os.chdir(folder+aml_file_name+"/")
    os.system("python loompy_to_readcount_indels.py -D {}".format(aml_file_name))
    seed=random.randint(1,random_number)
    with open("Seed.txt","w+") as f:
        f.write(str(seed))
    os.system("time python3 Algorithm_leiden_silhouette_on_indels.py -s "+ str(seed)+" -D "+ aml_file_name + "> output_log_" + aml_file_name + ".txt")
    #os.system("rm /home/priya/Downloads/Final_Stage/AML_Results/AML_remaining_datasets/"+aml_file_name+"/"+aml_file_name+".mpileup")
