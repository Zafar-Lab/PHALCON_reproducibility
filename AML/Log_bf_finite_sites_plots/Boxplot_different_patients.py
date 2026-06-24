import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# Read and merge batches
# =====================================================

df1 = pd.read_csv("Log_bf_exonic_batch1.tsv", sep="\t")
df2 = pd.read_csv("Log_bf_exonic_batch2.tsv", sep="\t")

combined_df = pd.concat([df1, df2], ignore_index=True)

combined_df.to_csv(
    "Log_bf_exonic_all_batches.tsv",
    sep="\t",
    index=False
)

# =====================================================
# Compute patient-level median log(BF)
# =====================================================

patient_summary = (
    combined_df
    .groupby("Data")["Log_BF_exonic_difference"]
    .median()
    .reset_index()
)

# Split into three equally sized groups
patient_summary["group"] = pd.qcut(
    patient_summary["Log_BF_exonic_difference"],
    q=3,
    labels=["Low", "Medium", "High"]
)

print("\nPatient group summary:")
print(
    patient_summary
    .groupby("group")["Log_BF_exonic_difference"]
    .agg(["count", "min", "median", "max"])
)

# =====================================================
# Add group labels back
# =====================================================

df = combined_df.merge(
    patient_summary[["Data", "group"]],
    on="Data"
)

# =====================================================
# Check number of patients in each group
# =====================================================

print("\nPatients per group:")
for grp in ["Low", "Medium", "High"]:
    n = df[df["group"] == grp]["Data"].nunique()
    print(grp, n)

# =====================================================
# Widths for final panels
# =====================================================

width_dict = {
    "Low": 5.5,
    "Medium": 5.5,
    "High": 5.5
}

# =====================================================
# Generate separate plots
# =====================================================

boxplot_color = "pink"
dot_color = "yellow"

for grp in ["Low", "Medium", "High"]:

    subset = df[df["group"] == grp].copy()

    # Sort AML_01_001, AML_02_001, ...
    subset["Data_sort_key"] = (
        subset["Data"]
        .str.extract(r"_(\d+)_")[0]
        .astype(int)
    )

    subset = subset.sort_values("Data_sort_key")

    sorted_order = subset["Data"].unique()

    plt.figure(figsize=(width_dict[grp], 10))

    sns.boxplot(
        data=subset,
        x="Data",
        y="Log_BF_exonic_difference",
        order=sorted_order,
        boxprops=dict(
            facecolor=boxplot_color,
            edgecolor="black",
            linewidth=0.8
        ),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
        medianprops=dict(color="black"),
        flierprops=dict(
            markerfacecolor="black",
            markeredgecolor="black"
        )
    )

    sns.stripplot(
        data=subset,
        x="Data",
        y="Log_BF_exonic_difference",
        order=sorted_order,
        color=dot_color,
        alpha=0.9,
        edgecolor="orange",
        linewidth=1,
        jitter=True
    )

   # plt.title(
   #     f"{grp} log(BF) patients "
   #     f"({subset['Data'].nunique()} patients)",
   #     fontsize=14
   # )

    plt.xlabel("", fontsize=12)
    plt.ylabel("", fontsize=12)

    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=16)

    plt.tight_layout()

    plt.savefig(
        f"Log_BF_{grp}_patients.png",
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()

print("\nFinished generating:")
print("Log_BF_Low_patients.png")
print("Log_BF_Medium_patients.png")
print("Log_BF_High_patients.png")