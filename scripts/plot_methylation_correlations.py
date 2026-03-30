#!/usr/bin/env python3
"""Plot a methylation correlation heatmap from pairwise results.

Usage: python plot_methylation_correlations.py <correlations_csv> <output_pdf>
Example: python plot_methylation_correlations.py results/all_pairs_correlations.csv figures/methylation_correlation_heatmap.pdf
"""
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

csv_path = sys.argv[1] if len(sys.argv) > 1 else sys.exit(f"Usage: {sys.argv[0]} <correlations_csv> <output_pdf>")
output_pdf = sys.argv[2]

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

MAPPING = {
    'EMSeq_HG002_LAB01_REP01.converted.bedGraph': 'HG002 WGBS',
    'EMSeq_HG003_LAB01_REP01.converted.bedGraph': 'HG003 WGBS',
    'EMSeq_HG004_LAB01_REP01.converted.bedGraph': 'HG004 WGBS',
    'CD_3033_GIAB.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG002 ONT Frozen',
    'CD_3029_GIAB.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG003 ONT Frozen',
    'CD_3031_GIAB.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG004 ONT Frozen',
    'CD_3032_Cache.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG002 ONT Ensilicated',
    'CD_3030_Cache.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG004 ONT Ensilicated',
    'CD_3028_Cache.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph': 'HG003 ONT Ensilicated',
}

ORDER = [
    'HG002 WGBS', 'HG003 WGBS', 'HG004 WGBS',
    'HG002 ONT Frozen', 'HG003 ONT Frozen', 'HG004 ONT Frozen',
    'HG002 ONT Ensilicated', 'HG003 ONT Ensilicated', 'HG004 ONT Ensilicated',
]

df = pd.read_csv(csv_path)
df['FileA_label'] = df['FileA'].map(MAPPING)
df['FileB_label'] = df['FileB'].map(MAPPING)

# build symmetric matrix from combinations (input may have each pair only once)
pivot_df = pd.DataFrame(np.nan, index=ORDER, columns=ORDER)
for _, row in df.iterrows():
    a, b, r = row['FileA_label'], row['FileB_label'], row['pearson_r']
    if a in ORDER and b in ORDER:
        pivot_df.loc[a, b] = r
        pivot_df.loc[b, a] = r
np.fill_diagonal(pivot_df.values, 1.0)

pivot_df.index.name = ''
pivot_df.columns.name = ''

mask = np.triu(np.ones_like(pivot_df, dtype=bool), k=1)

plt.figure(figsize=(7, 6))
sns.heatmap(pivot_df, mask=mask, annot=True, cmap='coolwarm',
            vmin=0.75, vmax=1, fmt=".2f")
plt.tight_layout()
plt.savefig(output_pdf)
print(f"Saved: {output_pdf}")
