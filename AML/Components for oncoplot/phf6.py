# Get all protein changes for PHF6
import pandas as pd
aml_variant_classification = pd.read_csv("variant_classification_aml_finite_sites.csv")
gene = "U2AF1"

phf6_df = (
    aml_variant_classification
    .query("Hugo_Symbol == @gene")
    [["Tumor_Sample_Barcode",
      "aaChange",
      "Variant_Classification",
      "Chromosome",
      "Start_Position"]]
    .sort_values("Tumor_Sample_Barcode")
)

print(phf6_df)