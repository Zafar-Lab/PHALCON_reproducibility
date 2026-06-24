import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# -------------------- Load data --------------------

phalcon = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/COMPASS vs GATK (tapestri insights) vs PHALCON vs Mutect2_all_samples AML_finite_sites/PHLACON_aml_precision_recall_on_mutect_ground_truth_finite_sites.tsv",
    sep='\t'
)
phalcon = phalcon[phalcon['AML name'].str.endswith('_001')]
print("number of 001 samples in phalcon:", phalcon.shape[0])

compass = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/COMPASS on GATK calls AML/compass_aml_precision_recall_on_mutect_ground_truth.tsv",
    sep='\t'
)
compass = compass[compass['AML name'].str.endswith('_001')]
print("number of 001 samples in compass:", compass.shape[0])

gatk = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/GATK_calls_AML/Gatk_aml_precision_recall_on_mutect_ground_truth.tsv",
    sep='\t'
)
gatk = gatk[gatk['AML name'].str.endswith('_001')]
print("number of 001 samples in gatk:", gatk.shape[0])

mutect2 = pd.read_csv(
    "Mutect2_aml_precision_recall_on_AML_ground_truth_all_samples.tsv",
    sep='\t'
)
mutect2 = mutect2[mutect2['AML_name'].str.endswith('_001')]
mutect2 = mutect2.rename(columns={'AML_name': 'AML name'})
print("number of 001 samples in mutect2:", mutect2.shape[0])

# -------------------- Add method labels --------------------
common_samples = set(phalcon['AML name']) & \
                 set(compass['AML name']) & \
                 set(gatk['AML name']) & \
                 set(mutect2['AML name'])

print("Number of common samples:", len(common_samples))

phalcon = phalcon[phalcon['AML name'].isin(common_samples)]
compass = compass[compass['AML name'].isin(common_samples)]
gatk = gatk[gatk['AML name'].isin(common_samples)]
mutect2 = mutect2[mutect2['AML name'].isin(common_samples)]

print("PHALCON:", phalcon.shape[0])
print("COMPASS:", compass.shape[0])
print("GATK:", gatk.shape[0])
print("Mutect2:", mutect2.shape[0])


phalcon['method'] = 'PHALCON'
compass['method'] = 'COMPASS'
gatk['method'] = 'GATK (Tapestri insights)'
mutect2['method'] = 'Mutect2'

# -------------------- Compute F1 --------------------

for df in [phalcon, compass, gatk, mutect2]:
    p = df['Precision']
    r = df['Recall']
    df['F1 score'] = np.where(
        (p + r) == 0,
        0,
        2 * p * r / (p + r)
    )

# -------------------- Combine --------------------

df_all = pd.concat([phalcon, compass, gatk, mutect2], ignore_index=True)
stats = (
    df_all.groupby('method')[['Precision','Recall','F1 score']]
    .agg(['mean','median','std'])
    .round(4)
)

print(stats)
# -------------------- Plot settings --------------------

method_order = ['COMPASS','GATK (Tapestri insights)','Mutect2','PHALCON']

xticks_annotation = [
    'COMPASS','GATK\n(Tapestri\ninsights)','Mutect2','PHALCON'
]

# ggplot-style theme
sns.set_theme(style="whitegrid", font_scale=1.2, rc={
    'axes.edgecolor': '0.8',
    'axes.linewidth': 1,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5
})

# SAME COLORS AS ORIGINAL SCRIPT
palette = {
    'COMPASS': '#EC8F96',
    'PHALCON': '#47CC69',
    'GATK (Tapestri insights)': '#EAD8A4',
    'Mutect2': "#77B1D4"
}

rc_parms = {
    "figure.figsize": [25, 6],
    "figure.dpi": 500,
    "font.size": 10,
    "font.family": "Arial"
}

save_parms = {"bbox_inches": "tight", "transparent": False}

with plt.rc_context(rc_parms):

    boxprops = dict(edgecolor='black', linewidth=1.)
    medianprops = dict(color='black', linewidth=2)
    whiskerprops = dict(color='black', linewidth=1.3)
    capprops = dict(color='black', linewidth=1.3)
    flierprops = dict(marker='o', markerfacecolor='black',
                      markeredgecolor='black', markersize=4, alpha=0.5)

    # -------------------- Plot grid --------------------

    fig, axes = plt.subplots(1, 3, figsize=(25, 6))

    for ax in axes:
        ax.set_facecolor('#e0e0e0')
        ax.tick_params(axis='both', colors='black')
        ax.yaxis.label.set_color('black')
        ax.xaxis.label.set_color('black')
        ax.title.set_color('black')
        ax.grid(color='white')

    metrics = ['Precision', 'Recall', 'F1 score']
    titles = ['Precision Comparison', 'Recall Comparison', 'F1 Score Comparison']

    for i, metric in enumerate(metrics):

        ax = axes[i]

        sns.boxplot(
            x='method',
            y=metric,
            data=df_all,
            ax=ax,
            width=0.4,
            order=method_order,
            palette=palette,
            boxprops=boxprops,
            medianprops=medianprops,
            whiskerprops=whiskerprops,
            capprops=capprops,
            flierprops=flierprops
        )

        ax.set_title(titles[i])
        ax.set_xlabel('')
        ax.set_ylabel(metric,fontsize=18)
        ax.set_xticklabels(xticks_annotation, rotation=0)
        ax.tick_params(axis='both', labelsize=15)

    # -------------------- Legend --------------------

    plt.subplots_adjust(right=0.82)

    legend_elements = [
        Patch(facecolor=palette['GATK (Tapestri insights)'], edgecolor='black', label='GATK (Tapestri insights)'),
        Patch(facecolor=palette['COMPASS'], edgecolor='black', label='COMPASS'),
        Patch(facecolor=palette['PHALCON'], edgecolor='black', label='PHALCON'),
        Patch(facecolor=palette['Mutect2'], edgecolor='black', label='Mutect2')
    ]

    fig.legend(
        handles=legend_elements,
        title='Method',
        loc='center right',
        bbox_to_anchor=(0.6, 1.01),
        frameon=True,
        fontsize=12,
        title_fontsize=13,
        ncol=4
    )

    plt.savefig(
        "PHALCON_finite_sites_vs_COMPASS_vs_GATK_vs_Mutect2_boxplot_AML_001_all_samples.png",
        **save_parms
    )
   