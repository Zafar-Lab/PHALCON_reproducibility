
# we first turn a gene's different mutations into different variants
# e.g. NPM1 has three variants 2234,2245,2267 then we break it up into NPM1_1, NPM1_2,NPM1_3
import os
import pandas as pd
import numpy as np
from collections import Counter

df_all_variants = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/Fraction of cells for manual oncoplot after srsf2 inclusion/AML_top_variants_fraction_of_mutated_cells.tsv",sep='\t')
# we remove intronic variants
df_all_variants = df_all_variants[df_all_variants['Intronic\Exonic']!='intronic']


samples = df_all_variants['Sample_no'].unique()
samples = samples.tolist() # list of all unique samples
samples=sorted(samples)
for unique_gene_name in set(df_all_variants['Gene Name']):  # extracting out the unique gene names
    print(unique_gene_name)
    df_for_that_gene = df_all_variants[df_all_variants['Gene Name'] == unique_gene_name]  # extracting the df of that particular gene
    frequency_of_each_variant = Counter(df_for_that_gene['Site'])  # Extracting out sites that are of that gene i.e. extracting all variants and finding their counts
    count = 1
    variants_in_decreasing_order_of_mutations = [key for key, count in frequency_of_each_variant.most_common()]
    print("variants in decreasing order of mutation:",variants_in_decreasing_order_of_mutations)
    for site in variants_in_decreasing_order_of_mutations: 
        df_all_variants.loc[df_all_variants['Site'] == site, 'Gene Name'] = unique_gene_name+"_"+str(count) # for each unique variant for that particular gene, we go to the original dataframe and update the gene name
        count+=1

df_all_variants.to_csv("Variant breakup fraction of cells across all AML samples.tsv",sep='\t',header=True)
df_all_variants = pd.read_csv("Variant breakup fraction of cells across all AML samples.tsv",sep='\t',index_col=[0])
variants_freq = Counter(df_all_variants['Gene Name'])
print("Frequency of variants :",variants_freq)
all_variants = [key for key, count in variants_freq.most_common()]  # extracting variant names e.g. NPM1_1, NPM1_2 and putting the into index, i.e. treating them as rows
manual_oncoplot = pd.DataFrame(0, index=all_variants, columns=samples) # columns are the samples

print("Making the dataframe for manual oncoplot ...")
for i in manual_oncoplot.index:  # variants names
    for j in manual_oncoplot.columns:  # sample names
        fraction = df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j), 'Fraction of cells mutated']
        #print(df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j)])
        
        try:
            manual_oncoplot.loc[i,j] = float(fraction) # if that variant is present in that sample, put the fraction
        except: 
            manual_oncoplot.loc[i,j] = 0 # if not, it will be an empty series and hence we put zero stating that this variant is not present in this sample


print(manual_oncoplot)
manual_oncoplot.to_csv("Dataframe_for_manual_oncoplot_creation_variant_breakup_AML.tsv",sep='\t',header=True,index=True)




df_all_variants = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/Fraction of cells for manual oncoplot after srsf2 inclusion/AML_top_variants_fraction_of_mutated_cells.tsv",sep='\t')



samples = df_all_variants['Sample_no'].unique()
samples = samples.tolist() # list of all unique samples
samples=sorted(samples)
for unique_gene_name in set(df_all_variants['Gene Name']):  # extracting out the unique gene names
    print(unique_gene_name)
    df_for_that_gene = df_all_variants[df_all_variants['Gene Name'] == unique_gene_name]  # extracting the df of that particular gene
    frequency_of_each_variant = Counter(df_for_that_gene['Site'])  # Extracting out sites that are of that gene i.e. extracting all variants and finding their counts
    count = 1
    variants_in_decreasing_order_of_mutations = [key for key, count in frequency_of_each_variant.most_common()]
    print("variants in decreasing order of mutation:",variants_in_decreasing_order_of_mutations)
    for site in variants_in_decreasing_order_of_mutations: # sites in decreasing order of the number of mutations in each gene
        variant_type = df_all_variants.loc[df_all_variants['Site'] == site, 'Intronic\Exonic']
        df_all_variants.loc[df_all_variants['Site'] == site, 'Gene Name'] = unique_gene_name+"_"+str(site)+"_"+variant_type # for each unique variant for that particular gene, we go to the original dataframe and update the gene name
        count+=1

df_all_variants.to_csv("Variant breakup fraction of cells across all AML samples site prefix.tsv",sep='\t',header=True)
df_all_variants = pd.read_csv("Variant breakup fraction of cells across all AML samples site prefix.tsv",sep='\t',index_col=[0])
variants_freq = Counter(df_all_variants['Gene Name'])
print("Frequency of variants :",variants_freq)
all_variants = [key for key, count in variants_freq.most_common()]  # extracting variant names e.g. NPM1_1, NPM1_2 and putting the into index, i.e. treating them as rows
manual_oncoplot = pd.DataFrame(0, index=all_variants, columns=samples) # columns are the samples

print("Making the dataframe for manual oncoplot ...")
for i in manual_oncoplot.index:  # variants names
    for j in manual_oncoplot.columns:  # sample names
        fraction = df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j), 'Fraction of cells mutated']
        #print(df_all_variants.loc[(df_all_variants['Gene Name'] == i) & (df_all_variants['Sample_no'] == j)])
        
        try:
            manual_oncoplot.loc[i,j] = float(fraction) # if that variant is present in that sample, put the fraction
        except: 
            manual_oncoplot.loc[i,j] = 0 # if not, it will be an empty series and hence we put zero stating that this variant is not present in this sample


print(manual_oncoplot)
manual_oncoplot.to_csv("Dataframe_for_manual_oncoplot_creation_variant_breakup_AML_site_prefix.tsv",sep='\t',header=True,index=True)

df = pd.read_csv("Dataframe_for_manual_oncoplot_creation_variant_breakup_AML_site_prefix.tsv",sep='\t',index_col=[0])
df[df > 0] = 1
no_of_samples = df.shape[1]
df['sum'] = df.sum(axis=1)
df.drop(df.columns[range(0,no_of_samples)], axis=1, inplace=True)
df.to_csv("Dataframe_for_manual_oncoplot_creation_AML_site_prefix_samples_affected_variant_breakup.tsv")



