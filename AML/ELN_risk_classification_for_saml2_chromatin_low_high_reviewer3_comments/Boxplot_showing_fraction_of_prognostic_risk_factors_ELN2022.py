import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================
# File paths
# =========================
file1 = "cluster1_mutations.csv"   # chromatin-low
file2 = "cluster2_mutations.csv"   # chromatin-high

# =========================
# Columns to sum
# =========================
risk_cols = [
    "ASXL1 presence", "BCOR presence", "EZH2 presence", "RUNX1 presence",
    "SF3B1 presence", "SRSF2 presence", "STAG2 presence", "U2AF1 presence",
    "ZRSR2 presence", "TP53 presence", "NPM1 presence and FLT3-ITD presence",
    "NPM1 absence and FLT3-ITD presence"
]

N_RISK_FACTORS = 12

# =========================
# Load and process
# =========================
def process_file(path, cluster_label):
    df = pd.read_csv(path, index_col=0)
    df["risk_sum"] = df[risk_cols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    df["risk_fraction"] = df["risk_sum"] 
    df["cluster"] = cluster_label
    print(df)
    return df[["risk_fraction", "cluster"]]

# =========================
# Combine clusters
# =========================
df_all = pd.concat([
    process_file(file1, "Chromatin mutated\ncell fraction low"),
    process_file(file2, "Chromatin mutated\ncell fraction high")
], ignore_index=True)

print(df_all)
# =========================
# Boxplot

sns.set(rc={'figure.figsize':(4,6)})
rc_parms = {"figure.figsize": [4, 6], "figure.dpi": 500, "font.size": 10, "font.family": "Arial"}
save_parms = {"bbox_inches": "tight", "transparent": True}
with plt.rc_context(rc_parms):
    sns.boxplot(data=df_all, x="cluster", y="risk_fraction", palette=["#CDEAC0", "#BFD7EA"])
    sns.stripplot(data=df_all, x="cluster", y="risk_fraction", color="black", alpha=0.6, jitter=0.15, size=4)

    # --- add mean as star ---
    means = df_all.groupby("cluster", sort=False)["risk_fraction"].mean()
    print(means)
    for i, mean_val in enumerate(means):
        plt.scatter(
            i,                  # x-position
            mean_val,           # y-position
            marker="*",         # star
            s=130,              # size
            color="#EFBF04",    # yellow
            edgecolor="black",
            linewidth=0.6,
            zorder=5
        )

    plt.ylabel("Number of ELN2022 risk factors present")
    plt.xticks(size=12)
    plt.tight_layout()
    plt.savefig('risk_factor_fraction_boxplot.png', bbox_inches='tight')

