#!/usr/bin/env python3
"""methylation calibration plots: observed EMSeq truth vs ONT ML probability.

bins single-read ML probabilities into deciles, computes observed methylation
fraction from EMSeq truth per bin, and plots calibration curves comparing
frozen and ensilicated conditions.

parallelized by chromosome, SSD-staged.

usage:
    python scripts/methylation_calibration.py

dependencies: pysam, numpy, pandas, matplotlib
"""
import os
import sys
import shutil
import numpy as np
import pandas as pd
import pysam
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

PROJECT = '/data/jb/project/giab_stanford'
METHYL_BAM_DIR = f'{PROJECT}/archive/data_pipeline/10_methylation_aligned'
MODKIT_DIR = f'{PROJECT}/archive/data_pipeline/11_modkit'
OUT_DIR = f'{PROJECT}/results/methylation_discordance'
FIG_DIR = f'{PROJECT}/figures'
SSD_STAGE = '/tmp/calib_stage'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(SSD_STAGE, exist_ok=True)

EMSEQ_METH_CUTOFF = 0.80
EMSEQ_UNMETH_CUTOFF = 0.20
REGIONS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
N_WORKERS = 23
N_BINS = 20  # number of probability bins

SAMPLES = {
    'HG002': {
        'frozen': 'CD_3033_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': 'CD_3032_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG002_LAB01_REP01.converted.bedGraph',
    },
    'HG003': {
        'frozen': 'CD_3031_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': 'CD_3030_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG003_LAB01_REP01.converted.bedGraph',
    },
    'HG004': {
        'frozen': 'CD_3029_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': 'CD_3028_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG004_LAB01_REP01.converted.bedGraph',
    },
}


def load_emseq_truth(emseq_file):
    """load EMSeq, return dict of chrom -> {pos -> truth_label}."""
    df = pd.read_csv(emseq_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'meth_pct', 'm_reads', 'u_reads'])
    df['coverage'] = df['m_reads'] + df['u_reads']
    df = df[(df['chrom'].isin(REGIONS)) & (df['coverage'] >= 10)]
    df['beta'] = df['meth_pct'] / 100.0

    truth = {}
    for chrom in REGIONS:
        cm = df[df['chrom'] == chrom]
        chrom_dict = {}
        for pos in cm[cm['beta'] >= EMSEQ_METH_CUTOFF]['start'].values:
            chrom_dict[int(pos)] = 1
        for pos in cm[cm['beta'] <= EMSEQ_UNMETH_CUTOFF]['start'].values:
            chrom_dict[int(pos)] = 0
        truth[chrom] = chrom_dict

    total = sum(len(v) for v in truth.values())
    print(f'  {os.path.basename(emseq_file)}: {total:,} sites', flush=True)
    return truth


def extract_chrom(args):
    """worker: extract (truth_label, ml_prob) pairs for one chromosome."""
    bam_path, chrom, chrom_truth = args
    labels = []
    probs = []

    bam = pysam.AlignmentFile(bam_path, 'rb')
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        try:
            mod_bases = read.modified_bases
        except (KeyError, ValueError):
            continue
        if not mod_bases:
            continue

        mC_key = None
        for key in mod_bases:
            if key[2] == 'm':
                mC_key = key
                break
        if mC_key is None:
            continue

        ref_positions = read.get_reference_positions(full_length=True)
        for query_pos, prob in mod_bases[mC_key]:
            if query_pos >= len(ref_positions):
                continue
            ref_pos = ref_positions[query_pos]
            if ref_pos is None:
                continue
            if ref_pos in chrom_truth:
                labels.append(chrom_truth[ref_pos])
                probs.append(prob)

    bam.close()
    return labels, probs


def extract_all(bam_path, truth):
    """extract all (truth, prob) pairs, parallelized by chromosome."""
    args = [(bam_path, chrom, truth.get(chrom, {})) for chrom in REGIONS]
    all_truth = []
    all_probs = []

    with Pool(N_WORKERS) as pool:
        for labels, probs in pool.imap_unordered(extract_chrom, args):
            all_truth.extend(labels)
            all_probs.extend(probs)

    return np.array(all_truth, dtype=np.int8), np.array(all_probs, dtype=np.uint8)


def compute_calibration(truth_arr, prob_arr, n_bins):
    """compute calibration curve: bin edges, observed fraction, predicted mean, counts."""
    prob_norm = prob_arr.astype(np.float64) / 255.0
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_idx = np.digitize(prob_norm, bin_edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    predicted = np.zeros(n_bins)
    observed = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=np.int64)

    for i in range(n_bins):
        mask = bin_idx == i
        n = mask.sum()
        counts[i] = n
        if n > 0:
            predicted[i] = prob_norm[mask].mean()
            observed[i] = truth_arr[mask].mean()

    return predicted, observed, counts


def stage_bam(bam_name):
    """copy BAM + CSI index to SSD, return SSD path."""
    hdd_bam = f'{METHYL_BAM_DIR}/{bam_name}'
    ssd_bam = f'{SSD_STAGE}/{bam_name}'
    shutil.copy2(hdd_bam, ssd_bam)
    for ext in ['.bai', '.csi']:
        hdd_idx = hdd_bam + ext
        if os.path.exists(hdd_idx):
            shutil.copy2(hdd_idx, ssd_bam + ext)
    return ssd_bam


def cleanup_bam(ssd_bam):
    """remove staged BAM + index."""
    os.remove(ssd_bam)
    for ext in ['.bai', '.csi']:
        idx = ssd_bam + ext
        if os.path.exists(idx):
            os.remove(idx)


def main():
    # collect calibration data per sample per condition
    all_calib = {}  # {sample: {cond: (predicted, observed, counts)}}
    all_tables = []

    for sample_name, paths in SAMPLES.items():
        print(f'\n=== {sample_name} ===', flush=True)
        truth = load_emseq_truth(paths['emseq'])
        all_calib[sample_name] = {}

        for cond in ['frozen', 'ensilicated']:
            bam_name = paths[cond]
            print(f'\n  staging {bam_name}...', flush=True)
            ssd_bam = stage_bam(bam_name)
            print(f'  extracting {cond}...', flush=True)

            truth_arr, prob_arr = extract_all(ssd_bam, truth)
            cleanup_bam(ssd_bam)

            print(f'  {len(truth_arr):,} calls, computing calibration...', flush=True)
            predicted, observed, counts = compute_calibration(truth_arr, prob_arr, N_BINS)
            all_calib[sample_name][cond] = (predicted, observed, counts)

            # save per-bin table
            for i in range(N_BINS):
                all_tables.append({
                    'sample': sample_name,
                    'condition': cond,
                    'bin': i,
                    'predicted_mean': round(predicted[i], 5),
                    'observed_fraction': round(observed[i], 5),
                    'n_calls': int(counts[i]),
                })

    # save table
    table_df = pd.DataFrame(all_tables)
    table_df.to_csv(f'{OUT_DIR}/calibration_data.tsv', sep='\t', index=False)

    # plot: one panel per sample, frozen + ensilicated overlaid
    colors = {'frozen': '#555555', 'ensilicated': '#00D8A4'}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    for ax, sample_name in zip(axes, SAMPLES.keys()):
        # perfect calibration line
        ax.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.4, label='perfect')

        for cond in ['frozen', 'ensilicated']:
            predicted, observed, counts = all_calib[sample_name][cond]
            mask = counts > 100  # skip empty bins
            label = cond.capitalize()
            ax.plot(predicted[mask], observed[mask], 'o-',
                    color=colors[cond], markersize=5, lw=1.5, label=label)

        ax.set_xlabel('Mean predicted probability')
        ax.set_ylabel('Observed fraction methylated')
        ax.set_title(sample_name, fontweight='bold')
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.legend(frameon=False, fontsize=9)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    fig.suptitle('Methylation calibration [EMSeq truth]', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/methylation_calibration.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/methylation_calibration.pdf', bbox_inches='tight')
    plt.close()

    # second figure: reliability diagram with histogram of predictions
    fig, axes = plt.subplots(2, 3, figsize=(13, 8),
                             gridspec_kw={'height_ratios': [3, 1]})

    for col, sample_name in enumerate(SAMPLES.keys()):
        ax_top = axes[0, col]
        ax_bot = axes[1, col]

        ax_top.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.4)

        for cond in ['frozen', 'ensilicated']:
            predicted, observed, counts = all_calib[sample_name][cond]
            mask = counts > 100
            label = cond.capitalize()
            ax_top.plot(predicted[mask], observed[mask], 'o-',
                        color=colors[cond], markersize=5, lw=1.5, label=label)

            # histogram of predictions
            bin_edges = np.linspace(0, 1, N_BINS + 1)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            ax_bot.bar(bin_centers, counts / counts.sum(),
                       width=1.0 / N_BINS * 0.4,
                       alpha=0.7, color=colors[cond],
                       align='center',
                       label=label if col == 0 else None)

        ax_top.set_xlim(-0.02, 1.02)
        ax_top.set_ylim(-0.02, 1.02)
        ax_top.set_aspect('equal')
        ax_top.set_title(sample_name, fontweight='bold')
        if col == 0:
            ax_top.set_ylabel('Observed fraction methylated')
        ax_top.legend(frameon=False, fontsize=9)
        for spine in ['top', 'right']:
            ax_top.spines[spine].set_visible(False)

        ax_bot.set_xlim(-0.02, 1.02)
        ax_bot.set_xlabel('Predicted probability')
        if col == 0:
            ax_bot.set_ylabel('Fraction of calls')
        for spine in ['top', 'right']:
            ax_bot.spines[spine].set_visible(False)

    fig.suptitle('Methylation reliability diagram [EMSeq truth]', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/methylation_reliability_diagram.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/methylation_reliability_diagram.pdf', bbox_inches='tight')
    plt.close()

    print(f'\nsaved to {FIG_DIR}/methylation_calibration.png')
    print(f'saved to {FIG_DIR}/methylation_reliability_diagram.png')
    print(f'saved to {OUT_DIR}/calibration_data.tsv')


if __name__ == '__main__':
    main()
