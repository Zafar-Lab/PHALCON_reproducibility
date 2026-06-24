import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

evolution_batch1 = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/Evolution_pattern_batch_1_finite_sites.txt",sep='\t')
evolution_batch2 = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/Evolution_pattern_batch_2_finite_sites.txt",sep='\t')

evolution_pattern = pd.concat([evolution_batch1,evolution_batch2])
print("Number of samples :",evolution_pattern.shape[0])
evolution_pattern.to_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/Evolution_pattern_barplot/Evolution_pattern_finite_sites_AML_exonic.txt",sep='\t',index=False)


evolution_pattern = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/Evolution_pattern_barplot/Evolution_pattern_finite_sites_AML_exonic.txt",sep='\t')


linear = evolution_pattern[evolution_pattern['Linear/Branching']=='Linear']
branching = evolution_pattern[evolution_pattern['Linear/Branching']=='Branching']

x = np.array(linear["Distinct genotypes"],dtype='int')
y = np.array(branching["Distinct genotypes"],dtype='int')


plt.hist([x,y],bins = np.arange(1,11) - 0.5, color = ['b','r'],label=['Linear','Branching'],alpha=0.4)
plt.legend()
plt.xlabel('Number of subclones')
plt.ylabel('Number of samples')
plt.title("Evolution pattern wrt number of distinct genotypes (PDX included) exonic")
plt.yticks(np.arange(0, 30.5, 5))
plt.xticks(range(1, 10))
#plt.savefig("Evolution pattern wrt number of distinct genotypes exonic with PDX finite sites")  no need as there are no pdx samples in finite sites already

plt.close()


# draw the same thing but remove the pdx samples
rows_to_drop = evolution_pattern['Sample_id'].str.contains('PDX')
from collections import Counter
evolution_pattern_without_pdx = evolution_pattern[~rows_to_drop]
evolution_pattern_without_pdx_001 = evolution_pattern_without_pdx[evolution_pattern_without_pdx['Sample_id'].str.endswith('_001')]
print("Evolution pattern without pdx 001 only : ",evolution_pattern_without_pdx_001)
print("Number of samples after pdx removal and 001 only : ",evolution_pattern_without_pdx_001.shape[0])
print("Count for linear vs branching :",Counter(evolution_pattern_without_pdx_001['Linear/Branching']))
linear = evolution_pattern_without_pdx_001[evolution_pattern_without_pdx_001['Linear/Branching']=='Linear']
branching = evolution_pattern_without_pdx_001[evolution_pattern_without_pdx_001['Linear/Branching']=='Branching']

x = np.array(linear["Distinct genotypes"],dtype='int')
y = np.array(branching["Distinct genotypes"],dtype='int')


plt.hist([x,y],bins = np.arange(1,11) - 0.5, color = ['b','r'],label=['Linear','Branching'],alpha=0.4)
plt.legend()
plt.xlabel('Number of subclones')
plt.ylabel('Number of samples')
plt.yticks(np.arange(0, 25.5, 5),fontsize=17)
plt.xticks(range(1, 10),fontsize=17)
plt.savefig("Evolution pattern wrt number of distinct genotypes exonic without PDX 001 samples only finite sites",dpi=300)