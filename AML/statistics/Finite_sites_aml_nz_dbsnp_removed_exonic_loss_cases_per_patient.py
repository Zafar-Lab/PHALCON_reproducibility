import os
import pandas as pd
frames=[]
# batch 1
BASE_DIR = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1"

all_exonic_loss = []
summary = []

# ---------------------------------------------------
# Helper: build genomic interval
# ---------------------------------------------------
def make_interval(df, chr_col, pos_col, ref_col, alt_col):
    df = df.copy()
    df["start"] = df[pos_col].astype(int)
    df["end"] = df["start"] + df[[ref_col, alt_col]].apply(
        lambda x: max(len(str(x[0])), len(str(x[1]))) - 1,
        axis=1
    )
    return df

# ---------------------------------------------------
# Iterate over patients
# ---------------------------------------------------
for patient in sorted(os.listdir(BASE_DIR)):
    patient_dir = os.path.join(BASE_DIR, patient)
    if not os.path.isdir(patient_dir):
        continue

    try:
        # -----------------------------
        # Locate files
        # -----------------------------
        geno_path = os.path.join(patient_dir, "Genotype configuration.tsv")

        final_df_file = next(
            f for f in os.listdir(patient_dir)
            if f.endswith("_indels_finite_sites_final_df.tsv")
        )
        var_func_file = next(
            f for f in os.listdir(patient_dir)
            if f.endswith("_indels_dbsnp_non_zero_finite_sites.avinput.variant_function")
        )

        final_df_path = os.path.join(patient_dir, final_df_file)
        var_func_path = os.path.join(patient_dir, var_func_file)

        # -----------------------------
        # Read files
        # -----------------------------
        geno = pd.read_csv(geno_path, sep="\t", header=None)
        final_df = pd.read_csv(final_df_path, sep="\t", header=None)

        var_func = pd.read_csv(
            var_func_path,
            sep="\t",
            header=None,
            names=[
                "region", "gene", "chr",
                "start", "end", "ref", "alt",
                "af", "dbsnp", "depth"
            ]
        )

        # -----------------------------
        # Sanity check
        # -----------------------------
        if geno.shape[0] != final_df.shape[0]:
            print(f"[SKIP] {patient}: row mismatch")
            continue

        # -----------------------------
        # Merge final_df + genotype
        # -----------------------------
        merged = final_df.iloc[:, :4].copy()
        merged.columns = ["chr", "site", "ref", "alt"]
        merged["config"] = geno.iloc[:, -1].values
        merged["patient"] = patient

        # -----------------------------
        # Keep only LOSS
        # -----------------------------
        loss_df = merged[merged["config"] == "Loss"].copy()

        if loss_df.empty:
            summary.append({
                "patient": patient,
                "num_exonic_loss": 0
            })
            print(f"[OK] {patient}: exonic loss = 0")
            continue

        # -----------------------------
        # Build intervals
        # -----------------------------
        loss_df = make_interval(loss_df, "chr", "site", "ref", "alt")

        exonic_var = var_func[var_func["region"] == "exonic"].copy()
        exonic_var = make_interval(exonic_var, "chr", "start", "ref", "alt")

        # -----------------------------
        # Interval overlap matching
        # -----------------------------
        matched_rows = []

        for _, l in loss_df.iterrows():
            hits = exonic_var[
                (exonic_var["chr"] == l["chr"]) &
                (exonic_var["start"] <= l["end"]) &
                (exonic_var["end"] >= l["start"])
            ]

            for _, h in hits.iterrows():
                matched_rows.append({
                    "chr": l["chr"],
                    "site": l["site"],
                    "ref": l["ref"],
                    "alt": l["alt"],
                    "gene": h["gene"],
                    "patient": patient
                })

        exonic_loss = pd.DataFrame(
            matched_rows,
            columns=["chr", "site", "ref", "alt", "gene", "patient"]
        )

        summary.append({
            "patient": patient,
            "num_exonic_loss": exonic_loss.shape[0]
        })

        if not exonic_loss.empty:
            all_exonic_loss.append(exonic_loss)

        print(f"[OK] {patient}: exonic loss = {exonic_loss.shape[0]}")

    except StopIteration:
        print(f"[SKIP] {patient}: required file missing")

    except Exception as e:
        print(f"[ERROR] {patient}: {e}")

# ---------------------------------------------------
# Final outputs
# ---------------------------------------------------
summary_df = pd.DataFrame(summary)

final_exonic_loss_df = (
    pd.concat(all_exonic_loss, ignore_index=True)
    if all_exonic_loss
    else pd.DataFrame(
        columns=["chr", "site", "ref", "alt", "gene", "patient"]
    )
)

summary_df.to_csv(
    os.path.join(BASE_DIR, "exonic_loss_summary.tsv"),
    sep="\t", index=False
)

final_exonic_loss_df.to_csv(
    os.path.join(BASE_DIR, "exonic_loss_variants_all_patients_AML_finite_sites.tsv"),
    sep="\t", index=False
)
frames.append(final_exonic_loss_df)
print("\nDONE.")


# batch 2
BASE_DIR = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2"

all_exonic_loss = []
summary = []

# ---------------------------------------------------
# Helper: build genomic interval
# ---------------------------------------------------
def make_interval(df, chr_col, pos_col, ref_col, alt_col):
    df = df.copy()
    df["start"] = df[pos_col].astype(int)
    df["end"] = df["start"] + df[[ref_col, alt_col]].apply(
        lambda x: max(len(str(x[0])), len(str(x[1]))) - 1,
        axis=1
    )
    return df

# ---------------------------------------------------
# Iterate over patients
# ---------------------------------------------------
for patient in sorted(os.listdir(BASE_DIR)):
    patient_dir = os.path.join(BASE_DIR, patient)
    if not os.path.isdir(patient_dir):
        continue

    try:
        # -----------------------------
        # Locate files
        # -----------------------------
        geno_path = os.path.join(patient_dir, "Genotype configuration.tsv")

        final_df_file = next(
            f for f in os.listdir(patient_dir)
            if f.endswith("_indels_finite_sites_final_df.tsv")
        )
        var_func_file = next(
            f for f in os.listdir(patient_dir)
            if f.endswith("_indels_dbsnp_non_zero_finite_sites.avinput.variant_function")
        )

        final_df_path = os.path.join(patient_dir, final_df_file)
        var_func_path = os.path.join(patient_dir, var_func_file)

        # -----------------------------
        # Read files
        # -----------------------------
        geno = pd.read_csv(geno_path, sep="\t", header=None)
        final_df = pd.read_csv(final_df_path, sep="\t", header=None)

        var_func = pd.read_csv(
            var_func_path,
            sep="\t",
            header=None,
            names=[
                "region", "gene", "chr",
                "start", "end", "ref", "alt",
                "af", "dbsnp", "depth"
            ]
        )

        # -----------------------------
        # Sanity check
        # -----------------------------
        if geno.shape[0] != final_df.shape[0]:
            print(f"[SKIP] {patient}: row mismatch")
            continue

        # -----------------------------
        # Merge final_df + genotype
        # -----------------------------
        merged = final_df.iloc[:, :4].copy()
        merged.columns = ["chr", "site", "ref", "alt"]
        merged["config"] = geno.iloc[:, -1].values
        merged["patient"] = patient

        # -----------------------------
        # Keep only LOSS
        # -----------------------------
        loss_df = merged[merged["config"] == "Loss"].copy()

        if loss_df.empty:
            summary.append({
                "patient": patient,
                "num_exonic_loss": 0
            })
            print(f"[OK] {patient}: exonic loss = 0")
            continue

        # -----------------------------
        # Build intervals
        # -----------------------------
        loss_df = make_interval(loss_df, "chr", "site", "ref", "alt")

        exonic_var = var_func[var_func["region"] == "exonic"].copy()
        exonic_var = make_interval(exonic_var, "chr", "start", "ref", "alt")

        # -----------------------------
        # Interval overlap matching
        # -----------------------------
        matched_rows = []

        for _, l in loss_df.iterrows():
            hits = exonic_var[
                (exonic_var["chr"] == l["chr"]) &
                (exonic_var["start"] <= l["end"]) &
                (exonic_var["end"] >= l["start"])
            ]

            for _, h in hits.iterrows():
                matched_rows.append({
                    "chr": l["chr"],
                    "site": l["site"],
                    "ref": l["ref"],
                    "alt": l["alt"],
                    "gene": h["gene"],
                    "patient": patient
                })

        exonic_loss = pd.DataFrame(
            matched_rows,
            columns=["chr", "site", "ref", "alt", "gene", "patient"]
        )

        summary.append({
            "patient": patient,
            "num_exonic_loss": exonic_loss.shape[0]
        })

        if not exonic_loss.empty:
            all_exonic_loss.append(exonic_loss)

        print(f"[OK] {patient}: exonic loss = {exonic_loss.shape[0]}")

    except StopIteration:
        print(f"[SKIP] {patient}: required file missing")

    except Exception as e:
        print(f"[ERROR] {patient}: {e}")

# ---------------------------------------------------
# Final outputs
# ---------------------------------------------------
summary_df = pd.DataFrame(summary)

final_exonic_loss_df = (
    pd.concat(all_exonic_loss, ignore_index=True)
    if all_exonic_loss
    else pd.DataFrame(
        columns=["chr", "site", "ref", "alt", "gene", "patient"]
    )
)

summary_df.to_csv(
    os.path.join(BASE_DIR, "exonic_loss_summary.tsv"),
    sep="\t", index=False
)

final_exonic_loss_df.to_csv(
    os.path.join(BASE_DIR, "exonic_loss_variants_all_patients_AML_finite_sites.tsv"),
    sep="\t", index=False
)
frames.append(final_exonic_loss_df)
print("\nDONE.")

pd.concat(frames).to_csv("all_exonic_loss_cases_aml_finite_sites_nz_dbsnp.tsv",sep='\t',index=False)
