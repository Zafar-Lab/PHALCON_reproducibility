import pandas as pd
import numpy as np

# -------------------- Load data --------------------

phalcon = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/COMPASS vs GATK (tapestri insights) vs PHALCON vs Mutect2_all_samples AML_finite_sites/PHLACON_aml_precision_recall_on_mutect_ground_truth_finite_sites.tsv",
    sep='\t'
)
phalcon = phalcon[phalcon['AML name'].str.endswith('_001')]

compass = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/COMPASS on GATK calls AML/compass_aml_precision_recall_on_mutect_ground_truth.tsv",
    sep='\t'
)
compass = compass[compass['AML name'].str.endswith('_001')]

gatk = pd.read_csv(
    "/home/priya/Downloads/Final_Stage/GATK_calls_AML/Gatk_aml_precision_recall_on_mutect_ground_truth.tsv",
    sep='\t'
)
gatk = gatk[gatk['AML name'].str.endswith('_001')]

mutect2 = pd.read_csv(
    "Mutect2_aml_precision_recall_on_AML_ground_truth_all_samples.tsv",
    sep='\t'
)
mutect2 = mutect2[mutect2['AML_name'].str.endswith('_001')]
mutect2 = mutect2.rename(columns={'AML_name': 'AML name'})

# -------------------- Keep only common samples --------------------

common_samples = (
    set(phalcon['AML name']) &
    set(compass['AML name']) &
    set(gatk['AML name']) &
    set(mutect2['AML name'])
)

phalcon = phalcon[phalcon['AML name'].isin(common_samples)].copy()
compass = compass[compass['AML name'].isin(common_samples)].copy()
gatk = gatk[gatk['AML name'].isin(common_samples)].copy()
mutect2 = mutect2[mutect2['AML name'].isin(common_samples)].copy()

# -------------------- Add method labels --------------------

phalcon['Method'] = 'PHALCON'
compass['Method'] = 'COMPASS'
gatk['Method'] = 'GATK'
mutect2['Method'] = 'Mutect2'

# -------------------- Function to compute F1 --------------------

def add_f1(df):
    p = df['Precision']
    r = df['Recall']

    df['F1'] = np.where(
        (p + r) == 0,
        0,
        2 * p * r / (p + r)
    )

    return df

phalcon = add_f1(phalcon)
compass = add_f1(compass)
gatk = add_f1(gatk)
mutect2 = add_f1(mutect2)

# -------------------- Select required columns --------------------

def format_df(df):
    return df[[
        'AML name',
        'Method',
        'Precision',
        'Recall',
        'F1'
    ]].rename(columns={
        'AML name': 'AML_sample_name'
    })

phalcon_final = format_df(phalcon)
compass_final = format_df(compass)
gatk_final = format_df(gatk)
mutect2_final = format_df(mutect2)

# -------------------- Combine all --------------------

final_df = pd.concat([
    phalcon_final,
    compass_final,
    gatk_final,
    mutect2_final
], ignore_index=True)

# -------------------- Sort nicely --------------------

final_df = final_df.sort_values(
    by=['AML_sample_name', 'Method']
).reset_index(drop=True)

# -------------------- Save --------------------

final_df.to_csv(
    "AML_precision_recall_f1_all_methods.tsv",
    sep='\t',
    index=False
)

###########################################################
# Pivot tables
###########################################################

pivot_precision = final_df.pivot(
    index='AML_sample_name',
    columns='Method',
    values='Precision'
)

pivot_recall = final_df.pivot(
    index='AML_sample_name',
    columns='Method',
    values='Recall'
)

pivot_f1 = final_df.pivot(
    index='AML_sample_name',
    columns='Method',
    values='F1'
)
###########################################################
# Function to compute percentage improvement
###########################################################

def percentage_improvement(df, metric_name):

    improvement_df = pd.DataFrame(index=df.index)

    for method in ['GATK', 'COMPASS', 'Mutect2']:

        baseline = df[method]

        # avoid division by zero
        improvement = np.where(
            baseline == 0,
            np.nan,
            ((df['PHALCON'] - baseline) / baseline) * 100
        )

        improvement_df[f'PHALCONfs_vs_{method}'] = improvement

    print(f"\n================ {metric_name} =================")

    avg_improvements = {}

    for col in improvement_df.columns:

        avg_val = np.nanmean(improvement_df[col])

        avg_improvements[col] = avg_val

        print(f"{col}: {avg_val:.2f}%")

    min_imp = min(avg_improvements.values())
    max_imp = max(avg_improvements.values())

    print(f"\nRange for {metric_name}: {min_imp:.2f}% - {max_imp:.2f}%")

    return improvement_df, avg_improvements
###########################################################
# Precision improvements
###########################################################

precision_improvement_df, precision_avg = percentage_improvement(
    pivot_precision,
    "Precision"
)

###########################################################
# Recall improvements
###########################################################

recall_improvement_df, recall_avg = percentage_improvement(
    pivot_recall,
    "Recall"
)

###########################################################
# F1 improvements
###########################################################

f1_improvement_df, f1_avg = percentage_improvement(
    pivot_f1,
    "F1 Score"
)