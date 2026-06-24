import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

chromatin_low = ['AML_02_001', 'AML_04_001', 'AML_07_001', 'AML_104_001', 'AML_105_001',
       'AML_118_001', 'AML_122_001', 'AML_39_001', 'AML_41_001', 'AML_55_001',
       'AML_75_001']
chromatin_high = ['AML_08_001', 'AML_103_001', 'AML_111_001', 'AML_114_001',
       'AML_116_001', 'AML_11_001', 'AML_18_001', 'AML_25_001', 'AML_29_001',
       'AML_42_001', 'AML_45_001', 'AML_46_001', 'AML_59_001', 'AML_60_001',
       'AML_67_001', 'AML_71_001', 'AML_77_001', 'AML_84_001', 'AML_87_001',
       'AML_89_001', 'AML_91_001', 'AML_94_001', 'AML_96_001', 'AML_98_001',
       'AML_99_001']

    
    
def get_asxl1(mut):
    if 'ASXL1' in mut:
        return 1
    else:
        return 0

def get_bcor(mut):
    if 'BCOR' in mut:
        return 1
    else:
        return 0
    
def get_ezh2(mut):
    if 'EZH2' in mut:
        return 1
    else:
        return 0

def get_runx1(mut):
    if 'RUNX1' in mut:
        return 1
    else:
        return 0

def get_sf3b1(mut):
    if 'SF3B1' in mut:
        return 1
    else:
        return 0

def get_srsf2(mut):
    if 'SRSF2' in mut:
        return 1
    else:
        return 0

def get_stag2(mut):
    if 'STAG2' in mut:
        return 1
    else:
        return 0

def get_u2af1(mut):
    if 'U2AF1' in mut or 'U2AF1;U2AF1L5' in mut:
        return 1
    else:
        return 0  

def get_zrsr2(mut):
    if 'ZRSR2' in mut:
        return 1
    else:
        return 0  

def get_tp53(mut):
    if 'TP53' in mut:
        return 1
    else:
        return 0  

def npm1_flt3(mut):
    if 'NPM1' in mut and 'FLT3-ITD' in mut: 
        return 1
    else:
        return 0  

def npm1_no_flt3(mut):
    if 'NPM1' not in mut and 'FLT3-ITD' in mut: 
        return 1
    else:
        return 0  


fractions=[]
frames=[]
genes=[]

for i in range(2):
    print(i)
    if i==0:
        cluster =chromatin_low
    if i==1:
        cluster=chromatin_high

    df=pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/ELN_risk_classification_for_saml2_chromatin_low_high_reviewer3_comments/Combined_AML_indels_nonzero_dbsnp_multianno_finite_sites.txt",sep='\t')
    
    df=df[~df['Func.refGene'].isin(['intronic','ncRNA_intronic','intergenic'])]
    # replace FLT3 with FLT3 itd wherever ref or alt is 0
    mask = (df['Ref'] == '0') & (df['Alt'] == '0')
    df['Gene.refGene'][mask] = 'FLT3-ITD'
 
    df=df[df['Tumor_Sample_Barcode'].isin(cluster)]
    df = df.groupby('Tumor_Sample_Barcode',as_index=False).agg(list)
    df = df[['Tumor_Sample_Barcode','Gene.refGene']]
    df['ASXL1 presence'] = df['Gene.refGene'].apply(get_asxl1)
    df['BCOR presence'] = df['Gene.refGene'].apply(get_bcor)
    df['EZH2 presence'] = df['Gene.refGene'].apply(get_ezh2)
    df['RUNX1 presence'] = df['Gene.refGene'].apply(get_runx1)
    df['SF3B1 presence'] = df['Gene.refGene'].apply(get_sf3b1)
    df['SRSF2 presence'] = df['Gene.refGene'].apply(get_srsf2)
    df['STAG2 presence'] = df['Gene.refGene'].apply(get_stag2)
    df['U2AF1 presence'] = df['Gene.refGene'].apply(get_u2af1)
    df['ZRSR2 presence'] = df['Gene.refGene'].apply(get_zrsr2)
    df['TP53 presence'] = df['Gene.refGene'].apply(get_tp53)
    df['NPM1 presence and FLT3-ITD presence'] = df['Gene.refGene'].apply(npm1_flt3)
    df['NPM1 absence and FLT3-ITD presence'] = df['Gene.refGene'].apply(npm1_no_flt3)

    cols=['ASXL1 presence','BCOR presence','EZH2 presence','RUNX1 presence','SF3B1 presence','SRSF2 presence','STAG2 presence','U2AF1 presence','ZRSR2 presence','TP53 presence','NPM1 presence and FLT3-ITD presence','NPM1 absence and FLT3-ITD presence']
    df['sum']=df[cols].sum(axis=1)
    df['sum'][df['sum']>0]=1
    genes_columns=cols
    frame = df[['Tumor_Sample_Barcode']+cols]
    frames.append(frame)
    
    one_cluster=[]
    for j in range(12):
        one_cluster.append(df[genes_columns[j]].sum())
    
    genes.append(one_cluster)

    print("at i:",i," genes:",genes[i])
    fraction=df['sum'].sum()/df.shape[0]
    fractions.append(fraction)
    print("Fraction of samples affected by at least one of these adverse prognostic risk factors for cluster"+str(i+1)+" is :"+str(round(fraction,4)))
    df.to_csv("cluster"+str(i+1)+"_mutations.csv")

colors=['#b598e4','#d798d3','#f29cc1','#f4bfb4','#f3e0a2']
colors_patch=["#d4bdb7",'#fd7f6f',"#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65","#beb9db","#fdcce5","#8bd3c7","#9381ff","#cb997e"]

genes=np.array(genes)
genes=genes.T
print(genes)
bottom = np.zeros(genes.shape[1])
for i in range(12):
        plt.bar(np.arange(2),genes[i],bottom=bottom,color=colors_patch[i])
        bottom += genes[i]
#colors_patch=["#A8D5BAFF",'#fd7f6f',"#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65","#beb9db","#fdcce5","#8bd3c7"]
patches=[]
idx = ['ASXL1 mut','BCOR mut','EZH2 mut','RUNX1 mut','SF3B1 mut','SRSF2 mut','STAG2 mut','U2AF1 mut','ZRSR2 mut','TP53 mut','NPM1 mut+FLT3-ITD mut','NPM1 wild type+FLT3-ITD mut']

for i in range(12):
    patches.append(mpatches.Patch(color=colors_patch[i], label=idx[i]))
plt.legend(handles=patches,bbox_to_anchor=(0, 1.35),loc=2,ncols=3)
plt.xticks(np.arange(2),['Chromatin low','Chromatin high'],fontsize=10)
plt.savefig("Fraction of samples affected by at least one of these adverse prognostic risk factors for each cluster ELN2022_adverse_intermediate",bbox_inches='tight')
plt.close()


df_concat = pd.concat(frames, axis=0)
print(df_concat)
df_concat.set_index('Tumor_Sample_Barcode',inplace=True)
df_concat =df_concat.transpose()
print(df_concat)


cluster_labels = [0]*len(chromatin_low)+[1]*len(chromatin_high)
custom_cluster_colors = {0:"#6D8764",1:"#E51400",2:"#1B9AAA",3:"#4B3F72",4:"#FFC43D"}
cluster_colors = [custom_cluster_colors[label] for label in cluster_labels]


# Convert custom colors to a valid colormap
cluster_cmap = ListedColormap([custom_cluster_colors[i] for i in range(5)])


no_color = '#edf2f4'
c_1= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#d4bdb7"])
c_2= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,'#fd7f6f'])
c_3= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#7eb0d5"])
c_4= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#b2e061"])
c_5= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#bd7ebe"])
c_6= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#ffb55a"])
c_7= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#ffee65"])
c_8= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#beb9db"])
c_9= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#fdcce5"])
c_10= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#8bd3c7"])
c_11= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#9381ff"])
c_12= matplotlib.colors.LinearSegmentedColormap.from_list("", [no_color,"#cb997e"])


cmap=sns.color_palette(cluster_colors, as_cmap=True)


cm = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10,c_11,c_12]
f, axs = plt.subplots(len(cm)+1, 1, figsize=(12, 4), gridspec_kw={'height_ratios': [0.3]+[1]*len(cm), 'hspace': 0})

sns.heatmap([cluster_labels], ax=axs[0], cmap=cluster_cmap, cbar=False)
axs[0].set_xticks([])  # Remove x-ticks
axs[0].set_yticks([])  # Remove y-ticks
axs[0].set_xticklabels([])  # Remove labels

counter = 0
for index, row in df_concat.iterrows():
    print(row)
    sns.heatmap(np.array([row.values]),  ax=axs[counter+1],cmap=cm[counter], cbar=False)
    #axs[counter+1].set_yticklabels([str(idx[counter])], rotation=0, ha='right')
    axs[counter+1].set_yticklabels([])
    axs[counter+1].set_xticklabels([])
    axs[counter+1].set_yticks([])
    counter += 1





plt.xticks([])
plt.yticks([])

#colors_patch=["#A8D5BAFF",'#fd7f6f',"#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65","#beb9db","#fdcce5","#8bd3c7"]
patches=[]
for i in range(12):
    patches.append(mpatches.Patch(color=colors_patch[i], label=idx[i]))


plt.legend(handles=patches,bbox_to_anchor=(-0.05, 15.5),loc=2,ncols=6)
plt.gcf().text(0.11,0.5, "Genetic abnormality for\n'Adverse' and 'Intermediate' risk category", ha="center", va="center", rotation=90, fontsize=10)
plt.xlabel("Samples")
plt.savefig("plot showing which of the four prognostic markers are active in which of the aml patients ELN2022_adverse_intermediate for chromatin low vs high patients",bbox_inches='tight')

print("chromatin low number of patients : ",len(chromatin_low))
print("chromatin high number of patients : ",len(chromatin_high))
