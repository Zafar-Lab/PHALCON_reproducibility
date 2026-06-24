import os
import pandas as pd

BASE_DIR = "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2/"

summary_rows = []
files=os.listdir(BASE_DIR)
for patient in files:
    pdir = os.path.join(BASE_DIR, patient)
    if not os.path.isdir(pdir):
        continue

    print(f"Processing {patient}")

    try:
        # ---------- final_df ----------
        final_df_file = [
            f for f in os.listdir(pdir)
            if f.endswith("_finite_sites_final_df.tsv")
        ][0]

        final_df = pd.read_csv(
            os.path.join(pdir, final_df_file),
            sep="\t",
            header=None
        ).iloc[:, :4]

        final_df.columns = ["chr", "site", "ref", "alt"]

        # ---------- Genotype configuration ----------
        geno = pd.read_csv(
            os.path.join(pdir, "Genotype configuration.tsv"),
            sep="\t",
            header=None
        )

        if len(final_df) != len(geno):
            print(f"[SKIP] {patient}: row mismatch")
            continue

        final_df["genotype"] = geno.iloc[:, -1].values
        loss_df = final_df[final_df["genotype"] == "Loss"]

        if loss_df.empty:
            continue

        # ---------- variant_function ----------
        vf_file = [
            f for f in os.listdir(pdir)
            if f.endswith(".avinput.variant_function")
        ][0]

        vf = pd.read_csv(
            os.path.join(pdir, vf_file),
            sep="\t",
            header=None
        )

        vf.columns = [
            "region", "gene", "chr", "start", "end",
            "ref", "alt"
        ] + list(vf.columns[7:])

        remove_regions = {"intronic", "ncRNA_intronic", "intergenic"}
        vf_exonic = vf[~vf["region"].isin(remove_regions)]

        matched_rows = []

        for _, lrow in loss_df.iterrows():
            chr_, pos, ref, alt = lrow["chr"], lrow["site"], lrow["ref"], lrow["alt"]

            # SNV or deletion
            if len(ref) >= len(alt):
                matches = vf_exonic[
                    (vf_exonic["chr"] == chr_) &
                    (vf_exonic["start"] == pos)
                ]

            # insertion
            else:
                matches = vf_exonic[
                    (vf_exonic["chr"] == chr_) &
                    (vf_exonic["start"] == pos) &
                    (alt[1:] == vf_exonic["alt"])
                ]
            
            #print(vf_exonic["start"] )
            #print(pos)

            if not matches.empty:
                for _, m in matches.iterrows():
                    matched_rows.append(m)
                    summary_rows.append(
                        [patient] + m.tolist()
                    )

        # ---------- Save per-patient ----------
        if matched_rows:
            out_df = pd.DataFrame(matched_rows, columns=vf.columns)
            out_df.to_csv(
                os.path.join(pdir, patient+"_LOSS_exonic_variants.tsv"),
                sep="\t",
                header=False,
                index=False
            )

            print(f"  → saved {len(out_df)} exonic loss variants")

    except Exception as e:
        print(f"[ERROR] {patient}: {e}")

# ---------- Global summary ----------
if summary_rows:
    summary_df = pd.DataFrame(
        summary_rows,
        columns=["patient"] + list(vf.columns)
    )

    summary_df.to_csv(
        os.path.join(BASE_DIR, "exonic_loss_summary.variant_function"),
        sep="\t",
        header=False,
        index=False
    )

print("DONE.")
