import os
import pandas as pd

frames = []
none = 0

def read_and_clean(path):
    df = pd.read_csv(path, sep="\t", header=None)

    # 🔴 DROP extra chr/pos columns if present
    # keep only first 9 columns (ANNOVAR format)
    if df.shape[1] > 9:
        df = df.iloc[:, :9]

    return df


# -------- Batch 1 --------
BASE_DIR = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/"
for file in os.listdir(BASE_DIR):
    try:
        variant_file = f"{BASE_DIR}/{file}/{file}_LOSS_exonic_variants.tsv"
        df = read_and_clean(variant_file)
        print(df)
        frames.append(df)
    except:
        none += 1
        print(file, "file not there")


# -------- Batch 2 --------
BASE_DIR = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/"
for file in os.listdir(BASE_DIR):
    try:
        variant_file = f"{BASE_DIR}/{file}/{file}_LOSS_exonic_variants.tsv"
        df = read_and_clean(variant_file)
        frames.append(df)
    except:
        none += 1
        print(file, "file not there")


# -------- Concatenate --------
concatenated_file = pd.concat(frames, ignore_index=True)

print("no loss case:", none)

concatenated_file.to_csv(
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/Final_summary_of_exonic_loss_variants.tsv",
    sep="\t",
    index=False,
    header=False
)
