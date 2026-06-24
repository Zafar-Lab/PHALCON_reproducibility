# intronic variants are removed here
import os
import pandas as pd
import numpy as np
from collections import Counter



df_all_variants = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/NPM1_and_saml2_further_stratification/AML_top_variants_fraction_of_mutated_cells_finite_sites.tsv",sep='\t')
for gene in ['FLT3','DNMT3A','NRAS','IDH2','NPM1','TET2','RUNX1','U2AF1','EZH2','PTPN11','BCOR','KRAS','WT1','ASXL1','SF3B1','IDH1','TP53','STAG2','PHF6','SETBP1','NF1','KIT','ZRSR2']:
    df_all_variants['Gene Name'] = df_all_variants['Gene Name'].apply(lambda x: gene if gene in x else x)

df_all_variants = df_all_variants[~df_all_variants['Intronic\Exonic'].str.contains('intronic')]
df_all_variants = df_all_variants[~df_all_variants['Intronic\Exonic'].str.contains('intergenic')]
df_all_variants = df_all_variants[~df_all_variants['Intronic\Exonic'].str.contains('ncRNA_intronic')]

print("Variants which are not intronic :\n",df_all_variants)
samples = df_all_variants['Sample_no'].unique()
samples = samples.tolist() # list of all unique samples
samples=sorted(samples)


variants_freq = Counter(df_all_variants['Gene Name'])
print("Frequency of variants :",variants_freq)
all_variants = [key for key, count in variants_freq.most_common()]  # extracting variant names e.g. NPM1_1, NPM1_2 and putting the into index, i.e. treating them as rows
manual_oncoplot = pd.DataFrame(0, index=all_variants, columns=samples) # columns are the samples

print("Making the dataframe for manual oncoplot ...")
for i in manual_oncoplot.index:  # variants names
    for j in manual_oncoplot.columns:  # sample names
        fraction = df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j), 'Fraction of cells mutated']
        #print(df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j)])
        if fraction.shape[0]>1:
            fraction=fraction.to_list()[0]
        try:
            manual_oncoplot.loc[i,j] = float(fraction) # if that variant is present in that sample, put the fraction
        except: 
            manual_oncoplot.loc[i,j] = 0 # if not, it will be an empty series and hence we put zero stating that this variant is not present in this sample




manual_oncoplot.to_csv("Pathway_Dataframe_for_manual_oncoplot_creation_AML_gene_wise_all_genes_from_all_pathways_included_finite_sites.tsv",sep='\t',header=True,index=True)


df = pd.read_csv("Pathway_Dataframe_for_manual_oncoplot_creation_AML_gene_wise_all_genes_from_all_pathways_included_finite_sites.tsv",sep='\t',index_col=[0])
df[df > 0] = 1
no_of_samples = df.shape[1]
print(df)
df['sum'] = df.sum(axis=1)
df = df.sort_values('sum',ascending=False)
df.drop(df.columns[range(0,no_of_samples)], axis=1, inplace=True)
df.to_csv("Gene_wise_AML_no_of_samples_affected.tsv")

manual_oncoplot = manual_oncoplot.reindex(df.index)
print(manual_oncoplot)
manual_oncoplot.to_csv("Pathway_Dataframe_for_manual_oncoplot_creation_AML_gene_wise_all_genes_from_all_pathways_included_finite_sites.tsv",sep='\t',header=True,index=True)

