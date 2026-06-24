import pandas as pd
import os

folder = '/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1/'

sub_folders = [
    name for name in os.listdir(folder)
    if os.path.isdir(os.path.join(folder, name))
]

frames = []

for aml_file_name in sub_folders:

    file_path = os.path.join(
        folder,
        aml_file_name,
        "Statistical_test_exonic",
        "Log BF test exonic.txt"
    )

    if os.path.exists(file_path):

        df = pd.read_csv(file_path, sep="\t")

        if {"Data",
            "Site number",
            "PHALCON_inferred_log_lklhd",
            "FP_inferred_log_lklhd"}.issubset(df.columns):

            df["Log_BF_exonic_difference"] = (
                df["PHALCON_inferred_log_lklhd"]
                - df["FP_inferred_log_lklhd"]
            )

            frames.append(
                df[["Data", "Site number", "Log_BF_exonic_difference"]]
            )

if frames:

    final_df = pd.concat(frames, ignore_index=True)
    log_bf_test = final_df[
    final_df["Data"].notna() &
    (final_df["Data"].astype(str) != "0.0")
]

    log_bf_test.to_csv(
        os.path.join(folder, "Log_bf_exonic_batch1.tsv"),
        sep="\t",
        index=False
    )

    print(f"Saved {len(log_bf_test)} rows")

else:
    print("No valid files found.")