import re
import os
import pandas as pd
import random
'''
This python file will go to each AML directory and run the files that are stated
'''

folder = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

frames = []

for aml_file_name in sub_folders:
    os.chdir(folder+aml_file_name+"/")
    evolution_pattern_file = pd.read_csv(aml_file_name+"_indels_evolution_pattern.txt",sep='\t',header=None)
    sample_name_evolution_pattern = evolution_pattern_file.values.tolist()[0]
    inferred_clusters_file = pd.read_csv("Final clusters inferred.txt",sep='\t')
    inferred_distinct_genotypes = inferred_clusters_file['Final inferred distinct genotypes'].values[0]

    sample_name_evolution_pattern.append(inferred_distinct_genotypes)

    frames.append(sample_name_evolution_pattern)
    

os.chdir(folder)
combined_evolution_pattern = pd.DataFrame(frames)
combined_evolution_pattern.rename(columns={0:"Sample_id",1:"Linear/Branching",2:"Distinct genotypes"},inplace=True)
combined_evolution_pattern.to_csv('Evolution_pattern_batch_1.txt',sep='\t',index=None) 