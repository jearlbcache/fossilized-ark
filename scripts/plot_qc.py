#!/usr/bin/env python3
"""Plot read length distributions comparing Frozen vs Ensilicated samples.

Usage: python plot_qc.py <fastq1> <fastq2> <label1> <label2> <output_prefix>
Example: python plot_qc.py archive/data_pipeline/13_concatenated_reads_batch2/CD_3745.fastq.gz \
             archive/data_pipeline/13_concatenated_reads_batch2/CD_3847.fastq.gz \
             "-80C" "Ensilicated" figures/3745_3847_batch2
"""
import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import SeqIO

mpl.rcParams.update({'font.family': 'Arial', 'font.size': 12})

fastq1 = sys.argv[1] if len(sys.argv) > 5 else sys.exit(
    f"Usage: {sys.argv[0]} <fastq1> <fastq2> <label1> <label2> <output_prefix>")
fastq2 = sys.argv[2]
label1 = sys.argv[3]
label2 = sys.argv[4]
output_prefix = sys.argv[5]


def read_lengths(fastq_path, max_records=int(1e7)):
    lengths = []
    with gzip.open(fastq_path, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if i >= max_records:
                break
            lengths.append(len(record.seq))
    return np.array(lengths)


def n50(lengths):
    s = np.sort(lengths)[::-1]
    return s[np.searchsorted(np.cumsum(s), s.sum() / 2)]


lengths1 = read_lengths(fastq1)
lengths2 = read_lengths(fastq2)

fig, ax = plt.subplots(figsize=(5, 5))
ax.hist(np.log10(lengths1), bins=50, density=True, alpha=1,
        label=f'{label1} (N50: {n50(lengths1)})', color='#DDE4E8')
ax.hist(np.log10(lengths2), bins=50, density=True, alpha=1,
        label=f'{label2} (N50: {n50(lengths2)})', color='#00D8A4')

for spine in ['top', 'left', 'right']:
    ax.spines[spine].set_visible(False)
ax.set_xlabel('Log10 read length')
ax.set_yticks([])
ax.legend(frameon=False, fontsize=10)

plt.tight_layout()
plt.savefig(f'{output_prefix}.png', dpi=300)
plt.savefig(f'{output_prefix}.pdf')
print(f"Saved: {output_prefix}.png/.pdf")
