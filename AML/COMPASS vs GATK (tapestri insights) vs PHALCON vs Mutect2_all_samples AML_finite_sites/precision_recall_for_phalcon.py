import pandas as pd
import numpy as np
import os


folder='/home/priya/Downloads/Final_Stage/AML_results_finite_sites/PHALCON_annovar_files/'
files=[i for i in os.listdir(folder)]
frames=[['AML name','Precision','Recall']]
for file in files: 
    aml_name = file[:-len('_indels_dbsnp_non_zero_finite_sites.avinput.variant_function')]
    print("aml name:",aml_name)
    print("File:",file)

    final_variants_df = pd.read_csv(folder+file,sep='\t',header=None)
    final_variants_df = final_variants_df[~final_variants_df[0].isin(['intronic','intergenic','ncRNA_intronic'])]
    final_variants_df=final_variants_df[[2,3,4,5,6]]
    final_variants_df.rename(columns={2:'chr',3:'pos',4:'pos_end',5:'ref',6:'alt'},inplace=True)
    print("phalcon variants df:\n",final_variants_df)

    ground_truth_df=pd.read_csv("/home/priya/Downloads/Final_Stage/GATK_calls_AML/Ground_truth_annovar/"+aml_name+"/"+aml_name+".avinput.variant_function",sep='\t',header=None)
    ground_truth_df = ground_truth_df[~ground_truth_df[0].isin(['intronic','intergenic','ncRNA_intronic'])]
    ground_truth_df=ground_truth_df[[2,3,4,5,6]]
    ground_truth_df.rename(columns={2:'chr',3:'pos',4:'pos_end',5:'ref',6:'alt'},inplace=True)
    print("ground truth df:\n",ground_truth_df)

    common = pd.merge(ground_truth_df, final_variants_df, on=['chr', 'pos', 'ref', 'alt'])
    print("Common variants:\n",common)

    TP=common.shape[0]
    FP=final_variants_df.shape[0]-common.shape[0]
    FN=ground_truth_df.shape[0]-common.shape[0]
    print("true positive:",TP)
    print("false positive:",FP)
    print("false negative:",FN)
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0


    frames.append([aml_name, precision, recall])

pd.DataFrame(frames).to_csv("PHLACON_aml_precision_recall_on_mutect_ground_truth_finite_sites.tsv",sep='\t',index=False,header=False)      