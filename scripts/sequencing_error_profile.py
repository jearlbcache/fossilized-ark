#!/usr/bin/env python3
"""per-read sequencing error profiles for frozen vs ensilicated BAMs.

extracts error rates from samtools stats and per-read identity from
CIGAR strings on chr20, comparing frozen and ensilicated conditions.

usage:
    python scripts/sequencing_error_profile.py

dependencies: pysam, numpy, pandas, matplotlib, subprocess
"""
import os
import subprocess
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

PROJECT = '/data/jb/project/giab_stanford'
BAM_DIR = f'{PROJECT}/archive/data_pipeline/12_giab_hg38_aligned'
OUT_DIR = f'{PROJECT}/results/error_profiles'
FIG_DIR = f'{PROJECT}/figures'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = {
    'HG002': {'frozen': 'CD_3033_GIAB', 'ensilicated': 'CD_3032_Cache'},
    'HG003': {'frozen': 'CD_3031_GIAB', 'ensilicated': 'CD_3030_Cache'},
    'HG004': {'frozen': 'CD_3029_GIAB', 'ensilicated': 'CD_3028_Cache'},
}

REGION = 'chr20'  # subsample to chr20 for per-read identity


def run_samtools_stats(bam_path):
    """run samtools stats and parse SN (summary numbers) lines."""
    print(f'    running samtools stats on {os.path.basename(bam_path)}...')
    result = subprocess.run(
        ['samtools', 'stats', bam_path],
        capture_output=True, text=True, env={**os.environ, 'PATH': '/home/anaconda3/bin:' + os.environ['PATH']})
    stats = {}
    for line in result.stdout.split('\n'):
        if line.startswith('SN\t'):
            parts = line.split('\t')
            key = parts[1].rstrip(':')
            try:
                val = float(parts[2])
            except (ValueError, IndexError):
                continue
            stats[key] = val
    return stats


def compute_per_read_identity(bam_path, region):
    """compute per-read identity from CIGAR on a region. returns array of identities."""
    print(f'    computing per-read identity on {region}...')
    bam = pysam.AlignmentFile(bam_path, 'rb')
    identities = []

    for read in bam.fetch(region):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        # from CIGAR: count matches (=), mismatches (X), insertions (I), deletions (D)
        # if CIGAR uses M (match/mismatch), use NM tag instead
        cigar = read.cigartuples
        if cigar is None:
            continue

        aligned_len = 0
        matches = 0
        has_eq = False

        for op, length in cigar:
            if op == 7:  # = (sequence match)
                matches += length
                aligned_len += length
                has_eq = True
            elif op == 8:  # X (sequence mismatch)
                aligned_len += length
                has_eq = True
            elif op == 0:  # M (match or mismatch)
                aligned_len += length
            elif op == 1:  # I (insertion)
                aligned_len += length
            elif op == 2:  # D (deletion)
                aligned_len += length

        if has_eq:
            # CIGAR uses =/X, we can compute identity directly
            if aligned_len > 0:
                identities.append(matches / aligned_len)
        else:
            # CIGAR uses M, use NM tag
            nm = read.get_tag('NM') if read.has_tag('NM') else None
            if nm is not None and aligned_len > 0:
                identities.append(1.0 - nm / aligned_len)

    bam.close()
    print(f'      {len(identities):,} reads')
    return np.array(identities)


def main():
    print('=== samtools stats (whole-genome) ===')
    error_results = []

    for sample_name, ids in SAMPLES.items():
        for cond in ['frozen', 'ensilicated']:
            bam_id = ids[cond]
            bam_path = f'{BAM_DIR}/{bam_id}.nist.aligned.sorted.bam'
            label = f'{sample_name}_{cond}'
            print(f'\n  {label}:')

            stats = run_samtools_stats(bam_path)

            bases_mapped = stats.get('bases mapped (cigar)', 0)
            mismatches = stats.get('mismatches', 0)
            insertions = stats.get('bases inserted', stats.get('insertions', 0))
            deletions = stats.get('bases deleted', stats.get('deletions', 0))

            mismatch_rate = mismatches / bases_mapped if bases_mapped > 0 else 0
            ins_rate = insertions / bases_mapped if bases_mapped > 0 else 0
            del_rate = deletions / bases_mapped if bases_mapped > 0 else 0

            print(f'    bases mapped: {bases_mapped:,.0f}')
            print(f'    mismatch rate: {mismatch_rate*100:.3f}%')
            print(f'    insertion rate: {ins_rate*100:.3f}%')
            print(f'    deletion rate: {del_rate*100:.3f}%')

            error_results.append({
                'sample': sample_name,
                'condition': cond,
                'bases_mapped': int(bases_mapped),
                'mismatches': int(mismatches),
                'insertions': int(insertions),
                'deletions': int(deletions),
                'mismatch_rate': round(mismatch_rate, 6),
                'insertion_rate': round(ins_rate, 6),
                'deletion_rate': round(del_rate, 6),
            })

    error_df = pd.DataFrame(error_results)
    error_df.to_csv(f'{OUT_DIR}/error_rates.tsv', sep='\t', index=False)
    print(f'\n{error_df.to_string(index=False)}')

    # per-read identity on chr20
    print(f'\n=== per-read identity ({REGION}) ===')
    identity_data = {}
    for sample_name, ids in SAMPLES.items():
        identity_data[sample_name] = {}
        for cond in ['frozen', 'ensilicated']:
            bam_path = f'{BAM_DIR}/{ids[cond]}.nist.aligned.sorted.bam'
            label = f'{sample_name}_{cond}'
            print(f'\n  {label}:')
            idents = compute_per_read_identity(bam_path, REGION)
            identity_data[sample_name][cond] = idents
            print(f'    median identity: {np.median(idents)*100:.2f}%')

    # plot per-read identity distributions
    colors = {'frozen': '#555555', 'ensilicated': '#00D8A4'}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    for ax, sample_name in zip(axes, SAMPLES.keys()):
        for cond in ['frozen', 'ensilicated']:
            idents = identity_data[sample_name][cond]
            med = np.median(idents) * 100
            ax.hist(idents * 100, bins=200, range=(80, 100), density=True,
                    histtype='step', lw=1.5, color=colors[cond],
                    label=f'{"Ensilicated" if cond == "ensilicated" else "Frozen"} [median {med:.1f}%]')
        ax.set_xlabel('Per-read identity [%]')
        if ax == axes[0]:
            ax.set_ylabel('Density')
        ax.set_title(sample_name, fontweight='bold')
        ax.legend(frameon=False, fontsize=9)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    fig.suptitle(f'Per-read identity distribution [{REGION}]', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/per_read_identity.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/per_read_identity.pdf', bbox_inches='tight')
    plt.close()
    print(f'\nsaved to {FIG_DIR}/per_read_identity.png')
    print(f'saved to {OUT_DIR}/error_rates.tsv')


if __name__ == '__main__':
    main()
