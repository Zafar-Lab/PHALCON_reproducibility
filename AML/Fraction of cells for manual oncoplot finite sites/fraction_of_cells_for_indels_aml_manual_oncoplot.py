# In order to create a manual oncoplot, we need to get the fraction of cells that are mutated 
# in the top genes as mentioned by the oncoplot generated using maftools

import os
import pandas as pd
import numpy as np

folder = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
#print(sub_folders)

print("Sub folders batch 1",sub_folders)

frames = []
frames.append(['Sample_no','Chromosome','Site','Ref','Alt','Fraction of cells mutated','Gene Name','Intronic\Exonic'])

for aml_sample_name in sub_folders:

   # try:
        vcf_obtained = np.loadtxt('/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/'+aml_sample_name+'/'+aml_sample_name+'_indels_nonzero_dbsnp_finite_sites_outputInference.vcf', dtype='str')
        print("vcf obtained :\n",vcf_obtained)
        print(vcf_obtained.shape[1]-9)
        vcf_df = pd.DataFrame(vcf_obtained)
        evidence_vcf = np.array(vcf_df)  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        print("VCF of evident sites :",evidence_vcf)
        annotated_file = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/'+ aml_sample_name + '/' +aml_sample_name + '_indels_dbsnp_non_zero_finite_sites.avinput.variant_function'
        annotated_df = pd.read_csv(annotated_file, sep = '\t', header = None)
        print("annotated df:",annotated_df)
        annotated_np_bulk = np.array(annotated_df,dtype='object')  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        
        genotype_arr = np.zeros((evidence_vcf.shape[0],evidence_vcf.shape[1]-9+2),dtype='object')


        for i in range(evidence_vcf.shape[0]):

            print(evidence_vcf[i])
            genotype_arr[i][0] = evidence_vcf[i][0]
            genotype_arr[i][1] = int(evidence_vcf[i][1])

            for j in range(evidence_vcf.shape[1]-9):
                if evidence_vcf[i][9:][j].split(":")[0] == '0/1':
                    genotype_arr[i][j+2] = 1
                elif evidence_vcf[i][9:][j].split(":")[0] == '0/0':
                    genotype_arr[i][j+2] = 0
                else:
                    print("Something wrong")

    
        print(genotype_arr)
        genotype_arr = pd.DataFrame(genotype_arr)
        
        print("Does genotype array of the evidence vcf have same shapes?:", genotype_arr.shape[0] == evidence_vcf.shape[0] == annotated_np_bulk.shape[0])
        print("genotype arrr:", genotype_arr)
        print("evidence vcf for bulk",evidence_vcf)

        numCells = genotype_arr.shape[1]-2
        cols  = list(range(2,numCells+2))
        genotype_arr['sum'] = genotype_arr[cols].sum(axis=1)

        print("new genotype arr\n",genotype_arr)

        for i in range(evidence_vcf.shape[0]):
                chrmsm = evidence_vcf[i][0]
                site = evidence_vcf[i][1]
                pos = evidence_vcf[i][1]
                ref = evidence_vcf[i][3]
                alt = evidence_vcf[i][4]
               
                gene_name = annotated_np_bulk[i][1]
                intronic_exonic = annotated_np_bulk[i][0]
                fraction_of_mutated_cells = genotype_arr.iloc[i]['sum']/numCells

                frames.append([aml_sample_name,chrmsm,site,ref,alt,fraction_of_mutated_cells, gene_name,intronic_exonic])
    
    #except:
     #   continue

folder = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/'

sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]
#print(sub_folders)

print("Sub folders batch 2",sub_folders)


for aml_sample_name in sub_folders:

   # try:
        vcf_obtained = np.loadtxt('/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/'+aml_sample_name+'/'+aml_sample_name+'_indels_nonzero_dbsnp_finite_sites_outputInference.vcf', dtype='str')
        print("vcf obtained :\n",vcf_obtained)
        print(vcf_obtained.shape[1]-9)
        vcf_df = pd.DataFrame(vcf_obtained)
        evidence_vcf = np.array(vcf_df)  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        print("VCF of evident sites :",evidence_vcf)
        annotated_file = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/'+ aml_sample_name + '/' +aml_sample_name + '_indels_dbsnp_non_zero_finite_sites.avinput.variant_function'
        annotated_df = pd.read_csv(annotated_file, sep = '\t', header = None)
        print("annotated df:",annotated_df)
        annotated_np_bulk = np.array(annotated_df,dtype='object')  # DONT THINK ABOUT THE VARIABLE NAME HERE, just used what was already there
        
        genotype_arr = np.zeros((evidence_vcf.shape[0],evidence_vcf.shape[1]-9+2),dtype='object')


        for i in range(evidence_vcf.shape[0]):

            print(evidence_vcf[i])
            genotype_arr[i][0] = evidence_vcf[i][0]
            genotype_arr[i][1] = int(evidence_vcf[i][1])

            for j in range(evidence_vcf.shape[1]-9):
                if evidence_vcf[i][9:][j].split(":")[0] == '0/1':
                    genotype_arr[i][j+2] = 1
                elif evidence_vcf[i][9:][j].split(":")[0] == '0/0':
                    genotype_arr[i][j+2] = 0
                else:
                    print("Something wrong")

    
        print(genotype_arr)
        genotype_arr = pd.DataFrame(genotype_arr)
        
        print("Does genotype array of the evidence vcf have same shapes?:", genotype_arr.shape[0] == evidence_vcf.shape[0] == annotated_np_bulk.shape[0])
        print("genotype arrr:", genotype_arr)
        print("evidence vcf for bulk",evidence_vcf)

        numCells = genotype_arr.shape[1]-2
        cols  = list(range(2,numCells+2))
        genotype_arr['sum'] = genotype_arr[cols].sum(axis=1)

        print("new genotype arr\n",genotype_arr)

        for i in range(evidence_vcf.shape[0]):
                chrmsm = evidence_vcf[i][0]
                site = evidence_vcf[i][1]
                pos = evidence_vcf[i][1]
                ref = evidence_vcf[i][3]
                alt = evidence_vcf[i][4]
               
                gene_name = annotated_np_bulk[i][1]
                intronic_exonic = annotated_np_bulk[i][0]
                fraction_of_mutated_cells = genotype_arr.iloc[i]['sum']/numCells

                frames.append([aml_sample_name,chrmsm,site,ref,alt,fraction_of_mutated_cells, gene_name,intronic_exonic])
    
    #except:
     #   continue
          

AML_top_variants_fraction_of_mutated_cells = pd.DataFrame(frames)
AML_top_variants_fraction_of_mutated_cells.to_csv('AML_top_variants_fraction_of_mutated_cells_finite_cells.tsv',sep='\t',index=None,header=False)    