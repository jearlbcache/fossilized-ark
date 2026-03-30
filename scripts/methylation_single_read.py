#!/usr/bin/env python3
"""single-read methylation concordance and ROC against EMSeq ground truth.

extracts per-read 5mC ML probabilities from ONT BAMs across all autosomes + chrX,
compares against EMSeq consensus truth at high-confidence sites
(beta >= 0.80 = methylated, beta <= 0.20 = unmethylated).
generates ROC curves, accuracy tables, and summary figures.

usage:
    python scripts/methylation_single_read.py

dependencies: pysam, numpy, pandas, sklearn, matplotlib
"""
import os
import sys
import numpy as np
import pandas as pd
import pysam
from sklearn.metrics import roc_curve, auc
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

METH_THRESHOLD = 128
EMSEQ_METH_CUTOFF = 0.80
EMSEQ_UNMETH_CUTOFF = 0.20

REGIONS = [f'chr{i}' for i in range(1, 23)] + ['chrX']

SAMPLES = {
    'HG002': {
        'frozen': f'{METHYL_BAM_DIR}/CD_3033_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': f'{METHYL_BAM_DIR}/CD_3032_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG002_LAB01_REP01.converted.bedGraph',
    },
    'HG003': {
        'frozen': f'{METHYL_BAM_DIR}/CD_3031_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': f'{METHYL_BAM_DIR}/CD_3030_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG003_LAB01_REP01.converted.bedGraph',
    },
    'HG004': {
        'frozen': f'{METHYL_BAM_DIR}/CD_3029_GIAB.methyl.aligned.sorted.bam',
        'ensilicated': f'{METHYL_BAM_DIR}/CD_3028_Cache.methyl.aligned.sorted.bam',
        'emseq': f'{MODKIT_DIR}/EMSeq_HG004_LAB01_REP01.converted.bedGraph',
    },
}


def load_emseq_truth(emseq_file):
    """load EMSeq genome-wide, return dict of chrom -> {pos -> truth_label}."""
    print(f'  loading EMSeq truth from {os.path.basename(emseq_file)}...')
    df = pd.read_csv(emseq_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'meth_pct', 'm_reads', 'u_reads'])
    df['coverage'] = df['m_reads'] + df['u_reads']
    df = df[(df['chrom'].isin(REGIONS)) & (df['coverage'] >= 10)]
    df['beta'] = df['meth_pct'] / 100.0

    truth = {}
    for chrom in REGIONS:
        cm = df[df['chrom'] == chrom]
        chrom_dict = {}
        meth = cm[cm['beta'] >= EMSEQ_METH_CUTOFF]
        unmeth = cm[cm['beta'] <= EMSEQ_UNMETH_CUTOFF]
        for pos in meth['start'].values:
            chrom_dict[int(pos)] = 1
        for pos in unmeth['start'].values:
            chrom_dict[int(pos)] = 0
        truth[chrom] = chrom_dict

    total = sum(len(v) for v in truth.values())
    total_meth = sum(sum(1 for x in v.values() if x == 1) for v in truth.values())
    print(f'    {total:,} high-confidence sites ({total_meth:,} meth, {total - total_meth:,} unmeth)')
    return truth


def extract_single_read_calls(bam_path, truth):
    """extract per-read 5mC calls at truth sites from a BAM, genome-wide."""
    print(f'  extracting reads from {os.path.basename(bam_path)}...')
    bam = pysam.AlignmentFile(bam_path, 'rb')

    total_truth = []
    total_probs = []
    reads_processed = 0

    for chrom in REGIONS:
        chrom_truth = truth.get(chrom, {})
        if not chrom_truth:
            continue

        chrom_labels = []
        chrom_probs = []

        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            reads_processed += 1

            try:
                mod_bases = read.modified_bases
            except (KeyError, ValueError):
                continue
            if not mod_bases:
                continue

            # find 5mC key
            mC_key = None
            for key in mod_bases:
                if key[2] == 'm':
                    mC_key = key
                    break
            if mC_key is None:
                continue

            # cache ref positions for this read (computed once)
            ref_positions = read.get_reference_positions(full_length=True)

            for query_pos, prob in mod_bases[mC_key]:
                if query_pos >= len(ref_positions):
                    continue
                ref_pos = ref_positions[query_pos]
                if ref_pos is None:
                    continue
                if ref_pos in chrom_truth:
                    chrom_labels.append(chrom_truth[ref_pos])
                    chrom_probs.append(prob)

        total_truth.extend(chrom_labels)
        total_probs.extend(chrom_probs)
        print(f'    {chrom}: {len(chrom_labels):,} calls ({reads_processed:,} reads cumulative)')
        sys.stdout.flush()

    bam.close()
    return np.array(total_truth, dtype=np.int8), np.array(total_probs, dtype=np.uint8)


def score_one_bam(label, bam_path, truth):
    """extract probabilities, compute accuracy and ROC, return result dict + ROC data."""
    print(f'\n  --- {label} ---')
    truth_arr, prob_arr = extract_single_read_calls(bam_path, truth)

    if len(truth_arr) == 0:
        print('  no calls extracted!')
        return None, None

    # accuracy at default threshold (128)
    call_arr = (prob_arr > METH_THRESHOLD).astype(np.int8)
    accuracy = (truth_arr == call_arr).mean()
    meth_mask = truth_arr == 1
    unmeth_mask = truth_arr == 0
    acc_meth = (call_arr[meth_mask] == 1).mean() if meth_mask.sum() > 0 else float('nan')
    acc_unmeth = (call_arr[unmeth_mask] == 0).mean() if unmeth_mask.sum() > 0 else float('nan')

    # ROC — subsample for sklearn (full dataset is too large for roc_curve)
    rng = np.random.default_rng(42)
    n_sub = min(5_000_000, len(truth_arr))
    idx = rng.choice(len(truth_arr), size=n_sub, replace=False)
    fpr, tpr, _ = roc_curve(truth_arr[idx], prob_arr[idx].astype(np.float32) / 255.0)
    roc_auc = auc(fpr, tpr)

    print(f'  overall accuracy:     {accuracy:.4f} ({len(truth_arr):,} calls)')
    print(f'  sensitivity (meth):   {acc_meth:.4f} ({meth_mask.sum():,} calls)')
    print(f'  specificity (unmeth): {acc_unmeth:.4f} ({unmeth_mask.sum():,} calls)')
    print(f'  AUC:                  {roc_auc:.4f}')

    result = {
        'condition': label,
        'region': 'genome-wide',
        'total_calls': len(truth_arr),
        'meth_calls': int(meth_mask.sum()),
        'unmeth_calls': int(unmeth_mask.sum()),
        'overall_accuracy': round(accuracy, 5),
        'sensitivity': round(acc_meth, 5),
        'specificity': round(acc_unmeth, 5),
        'auc': round(roc_auc, 5),
    }
    roc_data = {'fpr': fpr, 'tpr': tpr, 'auc': roc_auc}
    return result, roc_data


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    results = []
    # store ROC data per sample: {sample: {cond: roc_data}}
    all_rocs = {}

    for sample_name, paths in SAMPLES.items():
        print(f'\n=== {sample_name} ===')
        truth = load_emseq_truth(paths['emseq'])
        all_rocs[sample_name] = {}

        for cond in ['frozen', 'ensilicated']:
            label = f'{sample_name}_{cond}'
            r, roc_data = score_one_bam(label, paths[cond], truth)
            if r:
                results.append(r)
                all_rocs[sample_name][cond] = roc_data

    # save summary table
    summary = pd.DataFrame(results)
    summary.to_csv(f'{OUT_DIR}/single_read_concordance.tsv', sep='\t', index=False)
    print(f'\n{summary.to_string(index=False)}')

    # per-sample ROC plots (frozen and ensilicated overlaid)
    colors = {'frozen': '#555555', 'ensilicated': '#00D8A4'}
    cond_labels = {'frozen': 'Frozen', 'ensilicated': 'Ensilicated'}
    sample_names = list(SAMPLES.keys())

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, sample_name in zip(axes, sample_names):
        for cond in ['frozen', 'ensilicated']:
            rd = all_rocs[sample_name].get(cond)
            if rd is None:
                continue
            ax.plot(rd['fpr'], rd['tpr'], color=colors[cond], lw=1.5,
                    label=f'{cond_labels[cond]} [AUC={rd["auc"]:.4f}]')
        ax.plot([0, 1], [0, 1], 'k--', lw=0.5, alpha=0.3)
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(sample_name, fontweight='bold')
        ax.legend(loc='lower right', frameon=False, fontsize=9)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    fig.suptitle('Single-read methylation ROC [EMSeq truth]', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/methylation_single_read_roc.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/methylation_single_read_roc.pdf', bbox_inches='tight')
    plt.close()
    print(f'\nsaved ROC to {FIG_DIR}/methylation_single_read_roc.png')

    # zoomed ROC (top-left corner where the action is)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, sample_name in zip(axes, sample_names):
        for cond in ['frozen', 'ensilicated']:
            rd = all_rocs[sample_name].get(cond)
            if rd is None:
                continue
            ax.plot(rd['fpr'], rd['tpr'], color=colors[cond], lw=1.5,
                    label=f'{cond_labels[cond]} [AUC={rd["auc"]:.4f}]')
        ax.plot([0, 1], [0, 1], 'k--', lw=0.5, alpha=0.3)
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(sample_name, fontweight='bold')
        ax.legend(loc='lower right', frameon=False, fontsize=9)
        ax.set_xlim(-0.01, 0.25)
        ax.set_ylim(0.75, 1.01)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    fig.suptitle('Single-read methylation ROC [zoomed]', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/methylation_single_read_roc_zoomed.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/methylation_single_read_roc_zoomed.pdf', bbox_inches='tight')
    plt.close()
    print(f'saved zoomed ROC to {FIG_DIR}/methylation_single_read_roc_zoomed.png')

    print(f'\nsaved to {OUT_DIR}/single_read_concordance.tsv')


if __name__ == '__main__':
    main()
