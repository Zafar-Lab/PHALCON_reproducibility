# In order to create a manual oncoplot, we need to get the fraction of cells that are mutated 
# in the top genes as mentioned by the oncoplot generated using maftools
import pandas as pd
import numpy as np
import os

gene_groups = {
    "dna_meth": ['IDH1', 'IDH2', 'TET2', 'DNMT3A'],
    "npm1" : ['NPM1'],
    "rtk_ras":['NRAS','NF1','KRAS','FLT3','PTPN11','SETBP1','KIT','CBL','BRAF'],
    "tf":['RUNX1','WT1','PHF6','GATA2','ETV6'],
    "chromatin":['ASXL1','BCOR','EZH2','STAG2','KDM6A','SMC3','SMC1A'],
    "splicing":['SRSF2','U2AF1','SF3B1','ZRSR2'],
    "apoptosis" :['TP53','PPM1D','MYC','ATM','CHEK2','RAD21']
    }



def merge_redundant_rows(df):
    """
    Merge redundant rows based on gene name (column 4), perform OR operation on cell columns, 
    and keep the first row's information for the first 6 columns.
    """
    grouped = df.groupby(4, sort=False)

    merged_rows = []

    for gene, group in grouped:
        # Take the first row's non-cell columns
        non_cell_cols = group.iloc[0, :6].values

        # Perform OR operation across cell columns
        or_result = group.iloc[:, 6:].max().values

        # Combine non-cell columns with the OR results
        merged_row = list(non_cell_cols) + list(or_result)
        merged_rows.append((group.index[0], merged_row))  # Store the original index and row data

    # Sort merged rows by original index
    merged_rows.sort(key=lambda x: x[0])

    # Extract the rows (discarding indices) and create a new DataFrame
    merged_df = pd.DataFrame([row[1] for row in merged_rows], columns=df.columns)

    return merged_df


def group_genes_and_perform_or(df, gene_groups):
    """
    Perform OR operation within gene groups provided in a dictionary.

    Args:
        df (pd.DataFrame): The dataframe after merging duplicate rows.
        gene_groups (dict): Dictionary of gene groups.

    Returns:
        pd.DataFrame: Updated dataframe with grouped genes.
    """
    grouped_rows = []

    for group_name, genes in gene_groups.items():
        # Select rows belonging to the current group
        group_df = df[df[4].isin(genes)]

        if not group_df.empty:
            # Use the first row for non-cell columns
            non_cell_cols = group_df.iloc[0, :6].values

            # Perform OR operation across cell columns
            or_result = group_df.iloc[:, 6:].max().values

            # Combine non-cell columns with OR results
            grouped_row = list(non_cell_cols) + list(or_result)
            grouped_rows.append(grouped_row)

    # Append rows for genes that do not belong to any group
    ungrouped_genes = set(df[4]) - set(g for genes in gene_groups.values() for g in genes)
    ungrouped_df = df[df[4].isin(ungrouped_genes)]

    # Combine grouped and ungrouped rows
    final_df = pd.concat([pd.DataFrame(grouped_rows, columns=df.columns), ungrouped_df], ignore_index=True)

    return final_df


folder = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
#print(sub_folders)

print("Sub folders batch 1",sub_folders)

frames = []
frames.append(['Sample_no','Chromosome','Site','Ref','Alt','Fraction of cells mutated','Gene Name','Intronic\Exonic'])

for aml_sample_name in sub_folders:


        vcf_obtained = np.loadtxt(folder+aml_sample_name+'/'+aml_sample_name+'_indels_nonzero_dbsnp_finite_sites_outputInference.vcf', dtype='str')
        print("vcf obtained :\n",vcf_obtained)
        print(vcf_obtained.shape[1]-9)
        vcf_df = pd.DataFrame(vcf_obtained)
        evidence_vcf = np.array(vcf_df)  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        print("VCF of evident sites :",evidence_vcf)
        annotated_file = folder+ aml_sample_name + '/' +aml_sample_name + '_indels_dbsnp_non_zero_finite_sites.avinput.variant_function'
        annotated_df = pd.read_csv(annotated_file, sep = '\t', header = None)
        for gene in ['FLT3','DNMT3A','NRAS','IDH2','NPM1','TET2','RUNX1','U2AF1','EZH2','PTPN11','BCOR','KRAS','WT1','ASXL1','SF3B1','IDH1','TP53','STAG2','PHF6','SETBP1','NF1','KIT','ZRSR2']:
             annotated_df[1] = annotated_df[1].apply(lambda x: gene if gene in x else x)

        print("annotated df:",annotated_df)
        annotated_np_bulk = np.array(annotated_df,dtype='object')  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        
        genotype_arr = np.zeros((evidence_vcf.shape[0],evidence_vcf.shape[1]-9+6),dtype='object')


        for i in range(evidence_vcf.shape[0]):

            print(evidence_vcf[i])
            genotype_arr[i][0] = evidence_vcf[i][0]   # chrmsm
            genotype_arr[i][1] = int(evidence_vcf[i][1])  # site
            genotype_arr[i][2] = evidence_vcf[i][3]  # ref
            genotype_arr[i][3] = evidence_vcf[i][4]  # alt
            genotype_arr[i][4] = annotated_np_bulk[i][1]
            genotype_arr[i][5] = annotated_np_bulk[i][0] # intronic exonic

            for j in range(evidence_vcf.shape[1]-9):
                if evidence_vcf[i][9:][j].split(":")[0] == '0/1':
                    genotype_arr[i][j+6] = 1
                elif evidence_vcf[i][9:][j].split(":")[0] == '0/0':
                    genotype_arr[i][j+6] = 0
                else:
                    print("Something wrong")

    
        print(genotype_arr)
        genotype_arr = pd.DataFrame(genotype_arr)
        
        print("Does genotype array of the evidence vcf have same shapes?:", genotype_arr.shape[0] == evidence_vcf.shape[0] == annotated_np_bulk.shape[0])
        print("genotype arrr:", genotype_arr)
        print("evidence vcf for bulk",evidence_vcf)
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('ncRNA_intronic')]
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('intronic')]
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('intergenic')]
        genotype_arr = merge_redundant_rows(genotype_arr)
        genotype_arr = group_genes_and_perform_or(genotype_arr, gene_groups)


        numCells = genotype_arr.shape[1]-6
        cols  = list(range(6,numCells+6))
        genotype_arr['sum'] = genotype_arr[cols].sum(axis=1)
        

        print("new genotype arr\n",genotype_arr)
        


        for i in range(genotype_arr.shape[0]):
                chrmsm = genotype_arr.iloc[i][0]
                site = genotype_arr.iloc[i][1]
                ref = genotype_arr.iloc[i][2]
                alt = genotype_arr.iloc[i][3]
                gene_name = genotype_arr.iloc[i][4]
                intronic_exonic = genotype_arr.iloc[i][5]
                fraction_of_mutated_cells = genotype_arr.iloc[i]['sum']/numCells
                frames.append([aml_sample_name,chrmsm,site,ref,alt,fraction_of_mutated_cells, gene_name,intronic_exonic])

folder = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
#print(sub_folders)

print("Sub folders batch 2",sub_folders)


for aml_sample_name in sub_folders:


        vcf_obtained = np.loadtxt(folder+aml_sample_name+'/'+aml_sample_name+'_indels_nonzero_dbsnp_finite_sites_outputInference.vcf', dtype='str')
        print("vcf obtained :\n",vcf_obtained)
        print(vcf_obtained.shape[1]-9)
        vcf_df = pd.DataFrame(vcf_obtained)
        evidence_vcf = np.array(vcf_df)  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        print("VCF of evident sites :",evidence_vcf)
        annotated_file = folder+ aml_sample_name + '/' +aml_sample_name + '_indels_dbsnp_non_zero_finite_sites.avinput.variant_function'
        annotated_df = pd.read_csv(annotated_file, sep = '\t', header = None)
        for gene in ['FLT3','DNMT3A','NRAS','IDH2','NPM1','TET2','RUNX1','U2AF1','EZH2','PTPN11','BCOR','KRAS','WT1','ASXL1','SF3B1','IDH1','TP53','STAG2','PHF6','SETBP1','NF1','KIT','ZRSR2']:
            annotated_df[1] = annotated_df[1].apply(lambda x: gene if gene in x else x)

        print("annotated df:",annotated_df)
        annotated_np_bulk = np.array(annotated_df,dtype='object')  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        
        genotype_arr = np.zeros((evidence_vcf.shape[0],evidence_vcf.shape[1]-9+6),dtype='object')


        for i in range(evidence_vcf.shape[0]):

            print(evidence_vcf[i])
            genotype_arr[i][0] = evidence_vcf[i][0]   # chrmsm
            genotype_arr[i][1] = int(evidence_vcf[i][1])  # site
            genotype_arr[i][2] = evidence_vcf[i][3]  # ref
            genotype_arr[i][3] = evidence_vcf[i][4]  # alt
            genotype_arr[i][4] = annotated_np_bulk[i][1] # gene name
            genotype_arr[i][5] = annotated_np_bulk[i][0] # intronic exonic

            for j in range(evidence_vcf.shape[1]-9):
                if evidence_vcf[i][9:][j].split(":")[0] == '0/1':
                    genotype_arr[i][j+6] = 1
                elif evidence_vcf[i][9:][j].split(":")[0] == '0/0':
                    genotype_arr[i][j+6] = 0
                else:
                    print("Something wrong")

    
        print(genotype_arr)
        genotype_arr = pd.DataFrame(genotype_arr)
        
        print("Does genotype array of the evidence vcf have same shapes?:", genotype_arr.shape[0] == evidence_vcf.shape[0] == annotated_np_bulk.shape[0])
        print("genotype arrr:", genotype_arr)
        print("evidence vcf for bulk",evidence_vcf)
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('ncRNA_intronic')]
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('intronic')]
        genotype_arr = genotype_arr[~genotype_arr[5].str.contains('intergenic')]
        genotype_arr = merge_redundant_rows(genotype_arr)
        genotype_arr = group_genes_and_perform_or(genotype_arr, gene_groups)

        numCells = genotype_arr.shape[1]-6
        cols  = list(range(6,numCells+6))
        genotype_arr['sum'] = genotype_arr[cols].sum(axis=1)
        

        print("new genotype arr\n",genotype_arr)
        


        for i in range(genotype_arr.shape[0]):
                chrmsm = genotype_arr.iloc[i][0]
                site = genotype_arr.iloc[i][1]
                ref = genotype_arr.iloc[i][2]
                alt = genotype_arr.iloc[i][3]
                gene_name = genotype_arr.iloc[i][4]
                intronic_exonic = genotype_arr.iloc[i][5]
                fraction_of_mutated_cells = genotype_arr.iloc[i]['sum']/numCells
                frames.append([aml_sample_name,chrmsm,site,ref,alt,fraction_of_mutated_cells, gene_name,intronic_exonic])
    

          
header = frames[0]
AML_top_variants_fraction_of_mutated_cells = pd.DataFrame(frames[1:],columns=header)

# Replace gene names in column 4 with their corresponding group
gene_groups = {
    "dna_meth": ['IDH1', 'IDH2', 'TET2', 'DNMT3A'],
    "npm1" : ['NPM1'],
    "rtk_ras":['NRAS','NF1','KRAS','FLT3','PTPN11','SETBP1','KIT','CBL','BRAF'],
    "tf":['RUNX1','WT1','PHF6','GATA2','ETV6'],
    "chromatin":['ASXL1','BCOR','EZH2','STAG2','KDM6A','SMC3','SMC1A'],
    "splicing":['SRSF2','U2AF1','SF3B1','ZRSR2'],
    "apoptosis" :['TP53','PPM1D','MYC','ATM','CHEK2','RAD21']
    }

def invert_dict(d): 
    inverse = dict() 
    for key in d: 
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse: 
                # If not create a new list
                inverse[item] = key 
            else: 
                inverse[item].append(key) 
    return inverse

inverse= invert_dict(gene_groups)

AML_top_variants_fraction_of_mutated_cells["Gene Name"] = AML_top_variants_fraction_of_mutated_cells["Gene Name"].map(inverse).fillna(AML_top_variants_fraction_of_mutated_cells["Gene Name"])  # If a gene doesn't exist in the dictionary, keep it unchanged
print(AML_top_variants_fraction_of_mutated_cells)
AML_top_variants_fraction_of_mutated_cells.to_csv('AML_top_variants_fraction_of_mutated_cells_finite_sites.tsv',sep='\t',index=None,header=True)    
