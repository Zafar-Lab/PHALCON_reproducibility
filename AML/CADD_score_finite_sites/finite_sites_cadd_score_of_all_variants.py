import os
import pandas as pd

#########################################################
# PATHS
#########################################################

BIG_TABLE = "merged_cadd_aml_final_variants_srsf2_included_sample_wise_exonic_001_only.tsv"
VARFUNC_DIR = "finite_sites_variant_function_files_AML"
OUTPUT = "finite_sites_exonic_with_CADD.csv"

#########################################################
# READ BIG TABLE
#########################################################

big_df = pd.read_csv(BIG_TABLE,sep='\t')

# Rename if necessary
big_df = big_df.rename(columns={
    "sample no": "Sample_id",
    "site": "pos",
    "gene name": "gene",
    "cadd": "CADD"
})

big_df["chr"] = big_df["chr"].astype(str)
big_df["ref"] = big_df["ref"].astype(str)
big_df["alt"] = big_df["alt"].astype(str)

#########################################################
# EXACT LOOKUP
#########################################################

exact_lookup = {}

for idx, row in big_df.iterrows():

    key = (
        row["chr"],
        int(row["pos"]),
        row["ref"],
        row["alt"]
    )

    exact_lookup[key] = row["CADD"]

#########################################################
# ADD NEWLY COMPUTED CADD SCORES
#########################################################

extra_files = [
    "cadd_score_for_new_variants_set1.tsv",
    "cadd_score_for_new_variants_set2.tsv"
]

for f in extra_files:

    df = pd.read_csv(
        f,
        sep="\t",
        comment="#",
        header=None
    )

    # first 4 columns + last column(PHRED)
    df.columns = (
        ["Chrom", "Pos", "Ref", "Alt"]
        + [f"tmp{i}" for i in range(df.shape[1]-5)]
        + ["PHRED"]
    )

    for _, row in df.iterrows():

        key = (
            str(row["Chrom"]),
            int(row["Pos"]),
            str(row["Ref"]),
            str(row["Alt"])
        )

        exact_lookup[key] = row["PHRED"]

#########################################################
# FUNCTION TO FIND CADD
#########################################################

def find_cadd(chr_, pos, ref, alt):

    pos = int(pos)
    ref = str(ref)
    alt = str(alt)

    ##################################################
    # 1. Exact match
    ##################################################

    key = (chr_, pos, ref, alt)

    if key in exact_lookup:
        return exact_lookup[key]

    ##################################################
    # 2. Insertion:
    #
    #    variant_function:
    #
    #       - -> G
    #
    #    old:
    #
    #       C -> CG
    #
    ##################################################

    if ref in ["-", "*", "0"]:

        for _, r in big_df[
            (big_df["chr"] == chr_) &
            (big_df["pos"] == pos)
        ].iterrows():

            old_ref = r["ref"]
            old_alt = r["alt"]

            if (
                len(old_ref) == 1 and
                old_alt == old_ref + alt
            ):
                return r["CADD"]

    ##################################################
    # 3. Deletion:
    #
    # variant_function:
    #
    #    G -> -
    #
    # old:
    #
    #    AG -> A
    #
    ##################################################

    if alt in ["-", "*", "0"]:

        # same position
        for _, r in big_df[
            (big_df["chr"] == chr_) &
            (big_df["pos"] == pos)
        ].iterrows():

            old_ref = r["ref"]
            old_alt = r["alt"]

            if (
                len(old_alt) == 1 and
                old_ref == old_alt + ref
            ):
                return r["CADD"]

        # one base before
        for _, r in big_df[
            (big_df["chr"] == chr_) &
            (big_df["pos"] == pos - 1)
        ].iterrows():

            old_ref = r["ref"]
            old_alt = r["alt"]

            if (
                len(old_alt) == 1 and
                old_ref == old_alt + ref
            ):
                return r["CADD"]

    if ref == "0" and alt == "0":
        for (c, p, old_ref, old_alt), score in exact_lookup.items():
            if c == chr_ and p == pos:

                if old_ref.startswith("X") or old_alt.startswith("X"):
                    return score

        # otherwise if there is only one variant at this position
        vals = []

        for (c, p, _, _), score in exact_lookup.items():

            if c == chr_ and p == pos:
                vals.append(score)

        if len(vals) == 1:
            return vals[0]

    ##################################################
    # FINAL FALLBACK
    #
    # If there is exactly one variant at this
    # chromosome + position, use it.
    ##################################################

    vals = []

    for (c, p, _, _), score in exact_lookup.items():

        if c == chr_ and p == pos:
            vals.append(score)

    if len(vals) == 1:
        return vals[0]
    ##################################################
    # 4. Nothing found
    ##################################################

    return None


#########################################################
# PROCESS ALL VARIANT FUNCTION FILES
#########################################################

all_tables = []

for file in sorted(os.listdir(VARFUNC_DIR)):

    if not file.endswith(".variant_function"):
        continue

    sample = file.replace(
        "_indels_dbsnp_non_zero_finite_sites.avinput.variant_function",
        ""
    )

    path = os.path.join(VARFUNC_DIR, file)

    df = pd.read_csv(
        path,
        sep="\t",
        header=None
    )

    df.columns = [
        "annotation",
        "gene",
        "chr",
        "start",
        "end",
        "ref",
        "alt",
        "AF",
        "dummy1",
        "dummy2"
    ]

    ##################################################
    # Remove unwanted annotations
    ##################################################

    df = df[
        ~df["annotation"].isin(
            [
                "intronic",
                "intergenic",
                "ncRNA_intronic"
            ]
        )
    ].copy()

    ##################################################
    # Add Sample_id
    ##################################################

    df.insert(0, "Sample_id", sample)

    ##################################################
    # Find CADD
    ##################################################

    cadds = []

    for _, row in df.iterrows():

        cadd = find_cadd(
            row["chr"],
            row["start"],
            row["ref"],
            row["alt"]
        )

        cadds.append(cadd)

    df.insert(1, "CADD", cadds)

    all_tables.append(df)

#########################################################
# MERGE
#########################################################

final_df = pd.concat(
    all_tables,
    ignore_index=True
)

#########################################################
# SAVE
#########################################################

final_df.to_csv(
    OUTPUT,
    index=False
)



#########################################################
# REPORT
#########################################################
print("Total variants :", len(final_df))
print("Matched        :", final_df["CADD"].notna().sum())
print("Missing        :", final_df["CADD"].isna().sum())

print("\nMissing variants:\n")

print(
    final_df[
        final_df["CADD"].isna()
    ][
        [
            "Sample_id",
            "chr",
            "start",
            "ref",
            "alt",
            "gene"
        ]
    ]
)
