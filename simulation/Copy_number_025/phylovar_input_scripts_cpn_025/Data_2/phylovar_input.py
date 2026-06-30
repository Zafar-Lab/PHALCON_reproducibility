import pandas as pd
import numpy as np
data='2'
mpileup = pd.read_csv('sc_2000_output.mpileup',sep='\t',header=None)
print("Mpileup shape :",mpileup.shape)
vcf_obtained = np.loadtxt('/home/priya/Downloads/Final_Stage/Benchmarking/Copy_number_025/Sciphin_results/result'+data+'/data_cpn_025_'+data+'_2000_10.vcf',dtype='str')
filtered_sciphi_indices = np.array(vcf_obtained[:,1],dtype='int')
phylo_back_data_1_2000_10 = mpileup[mpileup[1].isin(filtered_sciphi_indices)]
print("Shape after filtering sciphi indices :", phylo_back_data_1_2000_10.shape)
phylo_back_data_1_2000_10.to_csv('phylo_input_data_cpn_025_'+data+'_2000_10.mpileup',sep='\t',header=False, index=False)