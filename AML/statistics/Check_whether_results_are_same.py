import os
import pandas as pd

folder_name = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_batch_1_finite_sites"

aml_samples = sorted([
    d for d in os.listdir(folder_name)
    if os.path.isdir(os.path.join(folder_name, d))
])

# counters
total = 0
skipped = 0
identical = 0
different = 0


for sample in aml_samples:
    finite_path   = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_batch_1_finite_sites/"+sample+"/Genotype configuration.tsv"
    infinite_path = "/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+sample+"/Genotype configuration.tsv"

    if not (os.path.exists(finite_path) and os.path.exists(infinite_path)):
        print(sample, "skip (missing files)")
        skipped += 1
        continue

    df1 = pd.read_csv(finite_path,sep='\t',header=None)
    df2 = pd.read_csv(infinite_path,sep='\t',header=None)

    last_row_equal = df1.iloc[-1].equals(df2.iloc[-1])
    print("Last row identical:", last_row_equal)

    

    # check labels + shape
    if not (df1.columns.equals(df2.columns) and df1.shape == df2.shape):
        print(sample, "skip (label/shape mismatch)")
        skipped += 1
        continue

    # compare values
    same = (df1 == df2).all().all()
    if same:
        print(sample, True)
        identical += 1
    else:
        print(sample, False)
        different += 1

# summary
print("\n===== SUMMARY =====")
print("Total samples:", total)
print("Skipped:", skipped)
print("Identical:", identical)
print("Different:", different)

for sample in aml_samples:
    df3 = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_batch_1_finite_sites/"+sample+"/"+sample+"_indels_finite_sites_final_df.tsv")
    df4 = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+sample+"/"+sample+"_indels_final_df.tsv")
    first4_equal = df3.iloc[:4].equals(df4.iloc[:4])
    print("First 4 rows identical:", first4_equal)
