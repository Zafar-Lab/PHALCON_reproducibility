import pandas as pd
import numpy as np
mpileup = pd.read_csv('sc_2000_output.mpileup',sep='\t',header=None)
vcf_obtained = np.loadtxt('data_2_2000_10.vcf',dtype='str')
filtered_sciphi_indices = np.array(vcf_obtained[:,1],dtype='int')
phylo_data_1_2000_10 = mpileup[mpileup[1].isin(filtered_sciphi_indices)]
phylo_data_1_2000_10.to_csv('phylo_input_data_2_2000_10.mpileup',sep='\t',header=False, index=False)