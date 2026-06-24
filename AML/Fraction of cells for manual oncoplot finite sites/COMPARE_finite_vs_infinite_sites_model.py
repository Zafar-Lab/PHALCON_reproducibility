import pandas as pd

inf = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/Fraction of cells for manual oncoplot after srsf2 inclusion/All non intronic variants (intronic, ncRNA_intronic, intergenic removed) non pdx samples aml", sep="\t")
fin = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/Fraction of cells for manual oncoplot finite sites/All non intronic variants (intronic, ncRNA_intronic, intergenic removed) non pdx samples aml finite sites.tsv", sep="\t")

# Ignore AML_002 samples
inf = inf[inf["Sample_no"].str.contains("001")]
fin = fin[fin["Sample_no"].str.contains("001")]


# Remove variants with zero fraction of cells mutated
inf = inf[inf["Fraction of cells mutated"] > 0]
fin = fin[fin["Fraction of cells mutated"] > 0]

count_inf = inf.groupby("Sample_no").size()
count_fin = fin.groupby("Sample_no").size()

comparison = pd.concat(
    [count_inf.rename("Infinite"),
     count_fin.rename("Finite")],
    axis=1
).fillna(0).astype(int)

comparison["Difference"] = comparison["Finite"] - comparison["Infinite"]

changed = comparison[comparison["Difference"] != 0]

print(changed)