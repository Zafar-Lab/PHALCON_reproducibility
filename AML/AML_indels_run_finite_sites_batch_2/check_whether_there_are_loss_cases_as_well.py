import pandas as pd
import os

base_dir = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/"

folders = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]
print(f"Total folders: {len(folders)}")

for folder in folders:
    folder_path = os.path.join(base_dir, folder)
    file_path = os.path.join(folder_path, "Genotype configuration.tsv")
    
    if not os.path.exists(file_path):
        print(f"⚠️ Missing file in {folder}")
        continue
    
    genotype_config = pd.read_csv(file_path, sep='\t', header=None)
    loss_geno = genotype_config[genotype_config.iloc[:, -1] == 'Loss']
    
    if not loss_geno.empty:
        print(f"\n==== {folder} ====")
        print(loss_geno)
