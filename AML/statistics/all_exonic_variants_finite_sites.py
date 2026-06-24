import os
import pandas as pd

BASE_DIRS = [
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1",
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2"
]

all_rows = []

for BASE_DIR in BASE_DIRS:

    # Iterate through each sample directory
    for sample_id in os.listdir(BASE_DIR):
        patient_dir = os.path.join(BASE_DIR, sample_id)

        if not os.path.isdir(patient_dir):
            continue

        try:
            # Find the variant function file
            var_func_file = next(
                f for f in os.listdir(patient_dir)
                if f.endswith("_indels_dbsnp_non_zero_finite_sites.avinput.variant_function")
            )
        except StopIteration:
            print(f"Variant function file not found for {sample_id}")
            continue

        var_func_path = os.path.join(patient_dir, var_func_file)

        # Read the variant function file
        df = pd.read_csv(var_func_path, sep="\t", header=None)

        # Keep only the required columns
        temp_df = pd.DataFrame({
            "Sample_id": sample_id,
            "chr": df[2],
            "start": df[3],
            "end": df[4],
            "ref": df[5],
            "alt": df[6],
            "gene_name": df[1].astype(str).str.split("(").str[0],
            "intronic/exonic": df[0]
        })

        # Remove unwanted categories
        temp_df = temp_df[
            ~temp_df["intronic/exonic"].isin(
                ["intronic", "ncRNA_intronic", "intergenic"]
            )
        ]

        all_rows.append(temp_df)

# Combine everything into a single dataframe
final_df = pd.concat(all_rows, ignore_index=True)

print(final_df)


# Optional: save
final_df.to_csv("AML_all_batches_variant_summary_finite_sites.csv", index=False)


flt3_zero = final_df[
    (final_df["gene_name"] == "FLT3") &
    (
        (final_df["ref"].astype(str) == "0") &
        (final_df["alt"].astype(str) == "0")
    )
]

print(flt3_zero)

flt3_zero_samples = sorted(flt3_zero["Sample_id"].unique())

print("Number of FLT3 samples with ref/alt = 0:", len(flt3_zero_samples))
print(flt3_zero_samples)