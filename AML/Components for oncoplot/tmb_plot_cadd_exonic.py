
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

background_color='#e2fdff'
bar_color='#f6de75'
color_of_points='#fc71c2'
sample_wise_variants_with_cadd_score = pd.read_csv('/home/priya/Downloads/Final_Stage/AML_results_finite_sites/CADD_score_finite_sites/finite_sites_exonic_with_CADD.csv')

tmb_df = sample_wise_variants_with_cadd_score.groupby('Sample_id').agg(
    total_entries=('Sample_id', 'size'),
    entries_greater_than_15=('CADD', lambda x: (x > 15).sum())
).reset_index()
tmb_df.columns=['sample no','number of variants','no. of variants with CADD>15']
print(sample_wise_variants_with_cadd_score)
print(tmb_df)
tmb_df.to_csv("cadd_tmb_for_oncoplot_finite_sites_number_of_samples.tsv",sep='\t',index=False)
