import pandas as pd
import numpy as np
mpileup = pd.read_csv('sc_2000_output.mpileup',sep='\t',header=None)
print("Mpileup shape :",mpileup.shape)
vcf_obtained = np.loadtxt('data_par_2_2000_10.vcf',dtype='str')
filtered_sciphi_indices = np.array(vcf_obtained[:,1],dtype='int')
phylo_par_data_2_2000_10 = mpileup[mpileup[1].isin(filtered_sciphi_indices)]
print("Shape after filtering sciphi indices :", phylo_par_data_2_2000_10.shape)
phylo_par_data_2_2000_10.to_csv('phylo_input_data_par_2_2000_10.mpileup',sep='\t',header=False, index=False)