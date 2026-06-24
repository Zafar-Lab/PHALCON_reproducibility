import os
import pandas as pd

# -----------------------------
# Change these paths
# -----------------------------
GROUND_TRUTH_DIR = "/home/priya/Downloads/Final_Stage/GATK_calls_AML/Ground_truth_annovar/"
PHALCON_BATCH1_folder = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/"
PHALCON_BATCH2_folder = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/"
OUTPUT_CSV = "finite_sites_exonic_not_detected_in_bulk.csv"

# Annotation types to ignore
REMOVE_TYPES = {
    "intronic",
    "intergenic",
    "ncRNA_intronic"
}

results = []
print(sorted(os.listdir(GROUND_TRUTH_DIR)))
# Iterate over sample folders
for sample in sorted(os.listdir(GROUND_TRUTH_DIR)):

    sample_folder = os.path.join(GROUND_TRUTH_DIR, sample)
    PHALCON_BATCH1 = os.path.join(PHALCON_BATCH1_folder, sample)
    PHALCON_BATCH2 = os.path.join(PHALCON_BATCH2_folder, sample)

    if not os.path.isdir(sample_folder):
        continue

    bulk_file = os.path.join(
        sample_folder,
        f"{sample}.avinput.variant_function"
    )
    
    
    phalcon_file = os.path.join(PHALCON_BATCH1, f"{sample}_indels_dbsnp_non_zero_finite_sites.avinput.variant_function")

    if not os.path.exists(phalcon_file):
        phalcon_file = os.path.join(PHALCON_BATCH2, f"{sample}_indels_dbsnp_non_zero_finite_sites.avinput.variant_function")
        print(phalcon_file)

    if not os.path.exists(phalcon_file):
        print(f"Missing PHALCON file for {sample}")
        continue
    
    if not os.path.exists(bulk_file):
        print(f"Missing bulk file for {sample}")
        continue


    # -----------------------------
    # Read files
    # -----------------------------
    bulk = pd.read_csv(
        bulk_file,
        sep="\t",
        header=None,
        dtype=str
    )

    phalcon = pd.read_csv(
        phalcon_file,
        sep="\t",
        header=None,
        dtype=str
    )

    # -----------------------------
    # Remove unwanted annotations
    # -----------------------------
    bulk = bulk[~bulk[0].isin(REMOVE_TYPES)].copy()
    phalcon = phalcon[~phalcon[0].isin(REMOVE_TYPES)].copy()

    # -----------------------------
    # Create comparison key
    # chr_start_end_ref_alt
    # columns:
    # 2 chr
    # 3 start
    # 4 end
    # 5 ref
    # 6 alt
    # -----------------------------
    bulk["key"] = (
        bulk[2] + "_" +
        bulk[3] + "_" +
        bulk[4] + "_" +
        bulk[5] + "_" +
        bulk[6]
    )

    phalcon["key"] = (
        phalcon[2] + "_" +
        phalcon[3] + "_" +
        phalcon[4] + "_" +
        phalcon[5] + "_" +
        phalcon[6]
    )

    bulk_keys = set(bulk["key"])

    unmatched = phalcon[~phalcon["key"].isin(bulk_keys)]

    unmatched_positions = unmatched[3].tolist()
    unmatched_positions = [int(i) for i in unmatched_positions]

    results.append({
        "Sample_id": sample,
        "PHALCON_inferred_sites": len(phalcon),
        "Sites that are not present in bulk": unmatched_positions
    })

# -----------------------------
# Save
# -----------------------------
df = pd.DataFrame(results)
df.to_csv(OUTPUT_CSV, index=False)

print(df)
print(f"\nSaved to {OUTPUT_CSV}")