import pandas as pd
from IPython.display import display
inf = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_Results/Fraction of cells for manual oncoplot after srsf2 inclusion/All non intronic variants (intronic, ncRNA_intronic, intergenic removed) non pdx samples aml", sep="\t")
fin = pd.read_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/Fraction of cells for manual oncoplot finite sites/All non intronic variants (intronic, ncRNA_intronic, intergenic removed) non pdx samples aml finite sites.tsv", sep="\t")

# Ignore AML_002 samples
inf = inf[inf["Sample_no"].str.contains("001")]
fin = fin[fin["Sample_no"].str.contains("001")]


# Gene sets per sample
inf_genes = inf.groupby("Sample_no")["Gene Name"].apply(set)
fin_genes = fin.groupby("Sample_no")["Gene Name"].apply(set)

# Find samples with different gene sets
changed_samples = []

for sample in sorted(set(inf["Sample_no"]).union(fin["Sample_no"])):
    genes_inf = inf_genes.get(sample, set())
    genes_fin = fin_genes.get(sample, set())

    if genes_inf != genes_fin:
        changed_samples.append(sample)

print(f"{len(changed_samples)} samples have different gene sets.")

# Print rows from both dataframes
for sample in changed_samples:
    print("\n" + "="*100)
    print(f"Sample: {sample}")

    genes_inf = inf_genes.get(sample, set())
    genes_fin = fin_genes.get(sample, set())

    print("Only in Infinite:", genes_inf - genes_fin)
    print("Only in Finite  :", genes_fin - genes_inf)

    print("\n--- Infinite ---")
    display(
        inf[inf["Sample_no"] == sample]
        .sort_values(["Gene Name", "Chromosome", "Site"])
    )

    print("\n--- Finite ---")
    display(
        fin[fin["Sample_no"] == sample]
        .sort_values(["Gene Name", "Chromosome", "Site"])
    )

# Gene sets per sample
inf_gene_sets = inf.groupby("Sample_no")["Gene Name"].apply(set)
fin_gene_sets = fin.groupby("Sample_no")["Gene Name"].apply(set)

all_samples = sorted(set(inf_gene_sets.index) | set(fin_gene_sets.index))

different_samples = []

for sample in all_samples:
    genes_inf = inf_gene_sets.get(sample, set())
    genes_fin = fin_gene_sets.get(sample, set())

    if genes_inf != genes_fin:
        different_samples.append(sample)

print(f"Number of samples with different gene sets: {len(different_samples)}")
print(different_samples)