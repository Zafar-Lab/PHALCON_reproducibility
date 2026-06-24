import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# reading the cadd dataframe for all variants (after srsf2 inclusion)
merged_df = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/CADD_score_finite_sites/finite_sites_exonic_with_CADD.csv")
#merged_df = merged_df.replace("Nan", 20)
final_variants_cadd_df = merged_df[merged_df['CADD'] != 'Nan'] # drop those variants where there is no cadd score

gene_order = ['FLT3', 'NRAS', 'DNMT3A', 'EZH2', 'NPM1', 'U2AF1', 'RUNX1', 'IDH2',
       'ASXL1', 'TET2', 'WT1', 'PTPN11', 'BCOR', 'KRAS', 'TP53', 'IDH1',
       'PHF6', 'SF3B1', 'SRSF2', 'KIT']
for gene in gene_order:
    final_variants_cadd_df['gene'] = final_variants_cadd_df['gene'].apply(lambda x: gene if gene in x else x)

final_variants_cadd_df.to_csv("finite_sites_exonic_with_CADD_top20.tsv",sep='\t',index=False)

##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
############################################################     PLOTTING TOP 20 variants            #########################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
##############################################################################################################################################################################################
df = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/CADD_score_finite_sites/finite_sites_exonic_with_CADD_top20.tsv",sep='\t')
frames = []

for gene in gene_order:
    gene_df = pd.DataFrame(df[df['gene'] == gene]['CADD'], index=None)
    gene_df.rename(columns={0:"CADD"},inplace=True)
    gene_df['Gene_name'] = gene
    #asxl1 = asxl1.drop_duplicates().reset_index().drop(columns='index').rename(columns={'cadd':'asxl1'})
    print(gene_df)
    frames.append(gene_df)

print(pd.concat(frames,axis=0))
pd.concat(frames,axis=0).to_csv('top_20_variants_boxplot_sample_wise_cadd_score_finite_sites_for_boxplot.tsv',sep='\t',index=False)

