#!/usr/bin/env python3
"""pairwise AUROC heatmap: every ONT BAM scored against every EMSeq truth.

scores all 6 ONT conditions (3 samples x frozen/ensilicated) against all 3
EMSeq truth sets, producing a 6x3 AUROC matrix rendered as a heatmap.

parallelized by chromosome within each BAM x truth pair.

usage:
    python scripts/methylation_pairwise_auroc.py

dependencies: pysam, numpy, pandas, sklearn, matplotlib, seaborn
"""
import os
import sys
import numpy as np
import pandas as pd
import pysam
from multiprocessing import Pool
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

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
SSD_STAGE = '/tmp/auroc_stage'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(SSD_STAGE, exist_ok=True)

EMSEQ_METH_CUTOFF = 0.80
EMSEQ_UNMETH_CUTOFF = 0.20
REGIONS = [f'chr{i}' for i in range(1, 23)] + ['chrX']
N_WORKERS = 23  # one per chromosome

# label -> (hdd_bam, bai)
BAM_FILES = {
    'HG002 frozen': 'CD_3033_GIAB.methyl.aligned.sorted.bam',
    'HG002 ensilicated': 'CD_3032_Cache.methyl.aligned.sorted.bam',
    'HG003 frozen': 'CD_3031_GIAB.methyl.aligned.sorted.bam',
    'HG003 ensilicated': 'CD_3030_Cache.methyl.aligned.sorted.bam',
    'HG004 frozen': 'CD_3029_GIAB.methyl.aligned.sorted.bam',
    'HG004 ensilicated': 'CD_3028_Cache.methyl.aligned.sorted.bam',
}

EMSEQ = {
    'HG002 EMSeq': f'{MODKIT_DIR}/EMSeq_HG002_LAB01_REP01.converted.bedGraph',
    'HG003 EMSeq': f'{MODKIT_DIR}/EMSeq_HG003_LAB01_REP01.converted.bedGraph',
    'HG004 EMSeq': f'{MODKIT_DIR}/EMSeq_HG004_LAB01_REP01.converted.bedGraph',
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
    total_meth = sum(sum(1 for x in v.values() if x == 1) for v in truth.values())
    print(f'  {os.path.basename(emseq_file)}: {total:,} sites ({total_meth:,} meth, {total - total_meth:,} unmeth)')
    return truth


def extract_chrom(args):
    """worker: extract single-read calls for one chromosome. returns (labels, probs)."""
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


def score_pair(bam_path, truth):
    """score one BAM against one truth set, parallelized by chromosome."""
    args = [(bam_path, chrom, truth.get(chrom, {})) for chrom in REGIONS]

    all_truth = []
    all_probs = []

    with Pool(N_WORKERS) as pool:
        for labels, probs in pool.imap_unordered(extract_chrom, args):
            all_truth.extend(labels)
            all_probs.extend(probs)

    if len(all_truth) == 0:
        return np.nan, 0

    truth_arr = np.array(all_truth, dtype=np.int8)
    prob_arr = np.array(all_probs, dtype=np.float32) / 255.0

    rng = np.random.default_rng(42)
    n_sub = min(5_000_000, len(truth_arr))
    idx = rng.choice(len(truth_arr), size=n_sub, replace=False)
    auc_val = roc_auc_score(truth_arr[idx], prob_arr[idx])

    return auc_val, len(truth_arr)


def main():
    print('loading EMSeq truth sets...')
    truths = {}
    for label, path in EMSEQ.items():
        truths[label] = load_emseq_truth(path)

    bam_labels = list(BAM_FILES.keys())
    truth_labels = list(EMSEQ.keys())
    results = []
    auc_matrix = np.full((len(bam_labels), len(truth_labels)), np.nan)

    for i, (bam_label, bam_name) in enumerate(BAM_FILES.items()):
        hdd_bam = f'{METHYL_BAM_DIR}/{bam_name}'
        ssd_bam = f'{SSD_STAGE}/{bam_name}'

        # stage BAM + index (CSI) to SSD
        import shutil
        print(f'\nstaging {bam_name} to SSD...', flush=True)
        shutil.copy2(hdd_bam, ssd_bam)
        for ext in ['.bai', '.csi']:
            hdd_idx = hdd_bam + ext
            if os.path.exists(hdd_idx):
                shutil.copy2(hdd_idx, ssd_bam + ext)
        print(f'  staged.', flush=True)

        for j, (truth_label, _) in enumerate(EMSEQ.items()):
            print(f'\n{bam_label} vs {truth_label}...', flush=True)
            auc_val, n_calls = score_pair(ssd_bam, truths[truth_label])
            auc_matrix[i, j] = auc_val
            results.append({
                'ont_condition': bam_label,
                'emseq_truth': truth_label,
                'auc': round(auc_val, 5),
                'n_calls': n_calls,
            })
            print(f'  AUC = {auc_val:.5f} ({n_calls:,} calls)', flush=True)

        # cleanup SSD
        os.remove(ssd_bam)
        for ext in ['.bai', '.csi']:
            idx = ssd_bam + ext
            if os.path.exists(idx):
                os.remove(idx)
        print(f'  cleaned up SSD.', flush=True)

    # save table
    res_df = pd.DataFrame(results)
    res_df.to_csv(f'{OUT_DIR}/pairwise_auroc.tsv', sep='\t', index=False)
    print(f'\n{res_df.to_string(index=False)}')

    # heatmap
    fig, ax = plt.subplots(figsize=(5, 6))
    annot = np.vectorize(lambda x: f'{x:.4f}')(auc_matrix)

    sns.heatmap(auc_matrix, annot=annot, fmt='',
                xticklabels=[l.replace(' EMSeq', '') for l in truth_labels],
                yticklabels=bam_labels,
                cmap='RdYlGn', vmin=0.80, vmax=0.95,
                linewidths=1.5, linecolor='white',
                cbar_kws={'label': 'AUROC'},
                ax=ax)
    ax.set_xlabel('EMSeq truth set')
    ax.set_ylabel('ONT condition')
    ax.set_title('Pairwise single-read methylation AUROC', fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/methylation_pairwise_auroc.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/methylation_pairwise_auroc.pdf', bbox_inches='tight')
    plt.close()
    print(f'\nsaved to {FIG_DIR}/methylation_pairwise_auroc.png')
    print(f'saved to {OUT_DIR}/pairwise_auroc.tsv')


if __name__ == '__main__':
    main()
