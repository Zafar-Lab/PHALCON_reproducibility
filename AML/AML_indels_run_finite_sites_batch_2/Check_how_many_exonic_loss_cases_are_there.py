import os
import pandas as pd

# CHANGE THIS to your main directory containing patient folders
base_dir = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2"

results = []

for patient in os.listdir(base_dir):
    patient_dir = os.path.join(base_dir, patient)

    if not os.path.isdir(patient_dir):
        continue

    try:
        # Find the avinput file
        avinput_files = [
            f for f in os.listdir(patient_dir)
            if f.endswith("_indels_non_zero_finite_sites.avinput")
        ]

        if len(avinput_files) == 0:
            print(f"[SKIP] {patient}: avinput file not found")
            continue

        avinput_path = os.path.join(patient_dir, avinput_files[0])
        geno_path = os.path.join(patient_dir, "Genotype configuration.tsv")

        # Read files
        av = pd.read_csv(avinput_path, sep="\t", header=None)
        geno = pd.read_csv(geno_path, sep="\t", header=None)

        # Sanity check: rows must match
        if av.shape[0] != geno.shape[0]:
            print(f"[SKIP] {patient}: row mismatch "
                  f"({av.shape[0]} vs {geno.shape[0]})")
            continue

        # Remove only intronic / ncRNA_intronic / intergenic
        remove_regions = ["intronic", "ncRNA_intronic", "intergenic"]
        keep_mask = ~av[0].isin(remove_regions)

        av_keep = av.loc[keep_mask].reset_index(drop=True)
        geno_keep = geno.loc[keep_mask].reset_index(drop=True)

        # Join row-wise
        joined = pd.concat([av_keep, geno_keep], axis=1)

        # Check for Loss in the last column
        genotype_col = joined.iloc[:, -1]
        loss_mask = genotype_col == "Loss"

        has_loss = loss_mask.any()
        num_loss = loss_mask.sum()

        results.append({
            "patient": patient,
            "has_loss": has_loss,
            "num_loss": num_loss
        })

        print(f"[OK] {patient}: Loss present = {has_loss}, count = {num_loss}")

    except Exception as e:
        print(f"[ERROR] {patient}: {e}")

# Create summary dataframe
summary_df = pd.DataFrame(results)

print("\n===== SUMMARY =====")
print(summary_df)

# Save summary
summary_df.to_csv("/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/loss_summary_all_patients.tsv", sep="\t", index=False)
