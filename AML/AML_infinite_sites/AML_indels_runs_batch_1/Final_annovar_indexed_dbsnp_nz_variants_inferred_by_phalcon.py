import re
import os
import pandas as pd

#####-------------EXPLANATION----------------####################
'''
We have variants from PHALCON (which have all zeros as well as dbsnp entries too)
We have final annovar output files which we obtain after removing all zeros and dbsnp
But annovar indexing is a little different in case of indels -> if "at position 1122 CAAG changes to C", annovar takes it as "at position 1123 AAG changes to * "
So, just so that we have a complete idea and a final readcount dataframe, i created this file which will make the final dataframe based on annovar indexing
'''

'''
This python file will go to each AML directory and make the dataframe of final variants
'''

folder = '/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/'
sub_folders = [name for name in os.listdir(folder) if os.path.isdir(os.path.join(folder, name))]

trues = []
for aml_file_name in sub_folders:
    print("File working on: ",aml_file_name)
    os.chdir(folder+aml_file_name+"/")
    phalcon_df = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+aml_file_name+"/"+aml_file_name+"_indels_final_df.tsv",sep='\t',header=None)
    final_variants = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/AML_indels_runs_batch_1/"+aml_file_name+"/"+aml_file_name+"_indels_dbsnp_non_zero.avinput.variant_function",header=None,sep='\t')
    phalcon_np = phalcon_df.to_numpy()
    for i in range(phalcon_np.shape[0]):
        if len(phalcon_np[i][2])>len(phalcon_np[i][3]) and phalcon_np[i][3]!="*":
           phalcon_np[i][1]+=1
           phalcon_np[i][2] = phalcon_np[i][2][1:]
           phalcon_np[i][3] = '*'
    phalcon_df = pd.DataFrame(phalcon_np) # annovar indexed dataframe
    dbsnp_nz_anno_phalcon_df = phalcon_df[phalcon_df[1].isin(final_variants[3])] # annovar indexed dataframe of final variants
    if dbsnp_nz_anno_phalcon_df.shape[0] == len(final_variants[3]): # is the shape of final variants matching to that reported by annovar
        trues.append(1)
    else:
        print("Problematic:",aml_file_name)

    dbsnp_nz_anno_phalcon_df.to_csv(aml_file_name+"_FINAL_VARIANTS_df_annovar_indexed_dbsnp_nz.tsv",sep='\t',header=False,index=False)

if sum(trues) == len(sub_folders):
    print("All files have variants equal to the annovar reported variants")
else:
    print("Problematic. Number of samples where it is not happening ",len(sub_folders)-sum(trues))