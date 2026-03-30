#!/usr/bin/env python3
"""Plot read length and Q-score distributions comparing two FASTQ files.

Usage: python plot_qscore_lengths.py <fastq1> <fastq2> <label1> <label2> <output_prefix>
Example: python plot_qscore_lengths.py archive/data_pipeline/03_concatenated_reads/CD_3033_GIAB.fastq \
             archive/data_pipeline/03_concatenated_reads/CD_3032_Cache.fastq \
             "-80C" "Ensilicated" figures/hg002_lengths_quality
"""
import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import SeqIO
from scipy.stats import ks_2samp

mpl.rcParams.update({'font.family': 'Arial', 'font.size': 12})

fastq1 = sys.argv[1] if len(sys.argv) > 5 else sys.exit(
    f"Usage: {sys.argv[0]} <fastq1> <fastq2> <label1> <label2> <output_prefix>")
fastq2 = sys.argv[2]
label1 = sys.argv[3]
label2 = sys.argv[4]
output_prefix = sys.argv[5]


def read_fastq(path, max_reads=None, qscore_subsample=1000):
    """Parse FASTQ, return (lengths, subsampled_qscores) arrays.

    Quality scores are subsampled (every Nth base) to avoid exhausting
    memory on large WGS datasets.
    """
    opener = gzip.open if path.endswith('.gz') else open
    lengths, qscores = [], []
    base_idx = 0
    with opener(path, "rt") as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if max_reads and i >= max_reads:
                break
            lengths.append(len(record.seq))
            for q in record.letter_annotations['phred_quality']:
                if base_idx % qscore_subsample == 0:
                    qscores.append(q)
                base_idx += 1
    return np.array(lengths), np.array(qscores)


def n50(lengths):
    s = np.sort(lengths)[::-1]
    return s[np.searchsorted(np.cumsum(s), s.sum() / 2)]


print(f"Reading {fastq1}...")
lengths1, qscores1 = read_fastq(fastq1)
print(f"Reading {fastq2}...")
lengths2, qscores2 = read_fastq(fastq2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

# read lengths
ax1.hist(np.log10(lengths1), bins=50, density=True, histtype='step',
         label=f'{label1} (N50: {n50(lengths1)})', color='k')
ax1.hist(np.log10(lengths2), bins=50, density=True, histtype='step',
         label=f'{label2} (N50: {n50(lengths2)})', color='#00D8A4')
ax1.set_xlabel('Log10 read length')
ax1.set_yticks([])
ax1.legend(frameon=False)
for spine in ['top', 'left', 'right']:
    ax1.spines[spine].set_visible(False)

# q-scores
ax2.hist(qscores1, bins=range(0, 41), density=True, histtype='step',
         label=label1, color='k')
ax2.hist(qscores2, bins=range(0, 41), density=True, histtype='step',
         label=label2, color='#00D8A4')
ax2.set_xlabel('Quality score [Phred]')
ax2.set_yticks([])
ax2.legend(frameon=False)
for spine in ['top', 'left', 'right']:
    ax2.spines[spine].set_visible(False)

ks_len, p_len = ks_2samp(lengths1, lengths2)
ks_q, p_q = ks_2samp(qscores1, qscores2)
print(f"K-S read lengths: statistic={ks_len:.4f}, p={p_len:.2e}")
print(f"K-S quality scores: statistic={ks_q:.4f}, p={p_q:.2e}")

plt.tight_layout()
plt.savefig(f'{output_prefix}.png', dpi=300)
plt.savefig(f'{output_prefix}.pdf')
print(f"Saved: {output_prefix}.png/.pdf")
