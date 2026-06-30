import pandas as pd
import numpy as np

data_no = 2
mpileup = pd.read_csv('sc_2000_output.mpileup',sep='\t',header=None)
print("Mpileup shape :",mpileup.shape)
vcf_obtained = np.loadtxt('data_fp_inc_'+str(data_no)+'_2000_10.vcf',dtype='str')
filtered_sciphi_indices = np.array(vcf_obtained[:,1],dtype='int')
phylo_data = mpileup[mpileup[1].isin(filtered_sciphi_indices)]
print("Shape after filtering sciphi indices :", phylo_data.shape)
phylo_data.to_csv('phylo_input_data_fp_inc_'+str(data_no)+'_2000_10.mpileup',sep='\t',header=False, index=False)