import pandas as pd
from collections import Counter

geno = pd.read_csv('Genotyped_dataframe_2000_50.csv',header=None)
print("Count of each case :",Counter(geno[2000]))