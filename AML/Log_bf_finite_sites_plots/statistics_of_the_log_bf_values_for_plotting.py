import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
df1 = pd.read_csv("Log_bf_exonic_batch1.tsv", sep="\t")
df2 = pd.read_csv("Log_bf_exonic_batch2.tsv", sep="\t")

combined_df = pd.concat([df1, df2], ignore_index=True)

combined_df.to_csv(
    "Log_bf_exonic_all_batches.tsv",sep="\t",index=False)

patient_summary = (
    combined_df.groupby('Data')['Log_BF_exonic_difference']
      .max()
      .reset_index()
      .sort_values('Log_BF_exonic_difference')
)

patient_summary.describe(
    percentiles=[0.25,0.5,0.75,0.9,0.95]
)

print(patient_summary)
