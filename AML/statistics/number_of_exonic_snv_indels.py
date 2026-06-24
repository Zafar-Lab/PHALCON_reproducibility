import os
import pandas as pd

base_dirs = [
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_1",
    "/home/priya/Downloads/Final_Stage/AML_results_finite_sites/AML_indels_run_finite_sites_batch_2"
]

allowed_annotations = {
    "exonic",
    "splicing",
    "UTR3",
    "UTR5"
}

rows = []

for base_dir in base_dirs:

    for sample in sorted(os.listdir(base_dir)):

        sample_dir = os.path.join(base_dir, sample)

        if not os.path.isdir(sample_dir):
            continue

        file_path = os.path.join(
            sample_dir,
            f"{sample}_indels_dbsnp_non_zero_finite_sites.avinput.variant_function"
        )

        if not os.path.exists(file_path):
            continue

        with open(file_path) as f:

            for line in f:

                cols = line.strip().split()

                if len(cols) < 7:
                    continue

                func = cols[0]

                if func not in allowed_annotations:
                    continue

                gene = cols[1]
                chrom = cols[2]
                start = cols[3]
                end = cols[4]
                ref = cols[5].upper()
                alt = cols[6].upper()


                rows.append([
    sample,
    func,
    gene,
    chrom,
    start,
    end,
    ref,
    alt
])

df = pd.DataFrame(
    rows,
    columns=[
        'Sample',
        'Func',
        'Gene',
        'Chr',
        'Start',
        'End',
        'Ref',
        'Alt'
    ]
)

df.to_csv('unique_variants_fs.tsv',sep='\t',index=False,header=False)
print("Before dedup:", len(df))

df_unique = df.drop_duplicates(
    subset=['Func','Gene','Chr','Start','End','Ref','Alt']
)
print("After dedup:", len(df_unique))
valid_bases = {'A','C','G','T'}

snv_mask = (
    df_unique['Ref'].isin(valid_bases) &
    df_unique['Alt'].isin(valid_bases)
)

indel_mask = ~snv_mask

print("Unique SNVs :", snv_mask.sum())
print("Unique Indels :", indel_mask.sum())

patients_with_flt3 = df.loc[
    df['Gene'] == 'FLT3',
    'Sample'
].nunique()

print("Patients with FLT3 mutation:", patients_with_flt3)
valid_bases = {'A','C','G','T'}

patients_with_flt3_indel = df.loc[
    (df['Gene'] == 'FLT3') &
    ~(
        df['Ref'].isin(valid_bases) &
        df['Alt'].isin(valid_bases)
    ),
    'Sample'
].nunique()
print("Patients with FLT3 indel/ITD:", patients_with_flt3_indel)
patients_with_no_flt3_indel = df.loc[
    (df['Gene'] == 'FLT3') &
    (
        df['Ref'].isin(valid_bases) &
        df['Alt'].isin(valid_bases)
    ),
    'Sample'
].nunique()

print("Patients with no FLT3 indel/ITD:", patients_with_no_flt3_indel)