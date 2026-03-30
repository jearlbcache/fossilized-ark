#!/usr/bin/env python3
"""discordant CpG methylation analysis between frozen and ensilicated samples.

compares per-CpG-site methylation between frozen and ensilicated (ensilicated)
conditions for GIAB Ashkenazi trio samples. identifies discordant sites
(|delta beta| > 0.20) and characterizes them by coverage, chromosome,
CpG island context, and GC content.

usage:
    python scripts/methylation_discordance.py

dependencies: pandas, numpy, matplotlib, pysam
"""
import os
import sys
import gzip
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

# configuration
PROJECT = '/data/jb/project/giab_stanford'
MODKIT_DIR = f'{PROJECT}/archive/data_pipeline/11_modkit'
REF_FASTA = f'{PROJECT}/reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta'
CPG_ISLAND_FILE = f'{PROJECT}/reference/cpgIslandExt.txt.gz'
OUT_DIR = f'{PROJECT}/results/methylation_discordance'
FIG_DIR = f'{PROJECT}/figures'

SAMPLES = {
    'HG002': {'frozen': 'CD_3033_GIAB', 'ensilicated': 'CD_3032_Cache'},
    'HG003': {'frozen': 'CD_3031_GIAB', 'ensilicated': 'CD_3030_Cache'},
    'HG004': {'frozen': 'CD_3029_GIAB', 'ensilicated': 'CD_3028_Cache'},
}

MIN_COV = 10
DISCORD_THRESH = 0.20
GC_WINDOW = 500
AUTOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX']

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


def load_bedgraph(sample_id):
    """load combined bedgraph: chrom, start, end, frac_modified, coverage."""
    path = f'{MODKIT_DIR}/{sample_id}.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph'
    print(f'  loading {os.path.basename(path)}...')
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'beta', 'coverage'],
                     dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                            'beta': np.float32, 'coverage': np.int32})
    return df


def load_cpg_islands():
    """load UCSC CpG island track, return intervals as arrays for fast lookup."""
    print('loading CpG island annotations...')
    islands = []
    with gzip.open(CPG_ISLAND_FILE, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom, start, end = parts[1], int(parts[2]), int(parts[3])
            islands.append((chrom, start, end))
    return pd.DataFrame(islands, columns=['chrom', 'start', 'end'])


def annotate_cpg_context(sites_df, islands_df):
    """annotate each site as island, shore [0-2kb], shelf [2-4kb], or open sea."""
    print('  annotating CpG context...')
    # build interval trees per chromosome for fast lookup
    from bisect import bisect_left, bisect_right

    island_by_chrom = defaultdict(list)
    for _, row in islands_df.iterrows():
        island_by_chrom[row['chrom']].append((row['start'], row['end']))

    # sort islands by start position
    for chrom in island_by_chrom:
        island_by_chrom[chrom].sort()

    labels = np.full(len(sites_df), 'open_sea', dtype=object)
    starts_arr = {ch: np.array([s for s, e in ivs]) for ch, ivs in island_by_chrom.items()}
    ends_arr = {ch: np.array([e for s, e in ivs]) for ch, ivs in island_by_chrom.items()}

    for chrom in sites_df['chrom'].unique():
        if chrom not in starts_arr:
            continue
        mask = sites_df['chrom'] == chrom
        positions = sites_df.loc[mask, 'start'].values
        isl_starts = starts_arr[chrom]
        isl_ends = ends_arr[chrom]
        chrom_labels = np.full(len(positions), 'open_sea', dtype=object)

        for i, pos in enumerate(positions):
            # find nearest island
            idx = bisect_right(isl_starts, pos) - 1
            min_dist = float('inf')

            for j in [idx - 1, idx, idx + 1]:
                if 0 <= j < len(isl_starts):
                    if isl_starts[j] <= pos < isl_ends[j]:
                        min_dist = 0
                        break
                    dist = min(abs(pos - isl_starts[j]), abs(pos - isl_ends[j]))
                    min_dist = min(min_dist, dist)

            if min_dist == 0:
                chrom_labels[i] = 'island'
            elif min_dist <= 2000:
                chrom_labels[i] = 'shore'
            elif min_dist <= 4000:
                chrom_labels[i] = 'shelf'

        labels[mask.values] = chrom_labels

    return labels


def compute_gc_content(sites_df, ref):
    """compute GC fraction in a window centered on each site."""
    print('  computing GC content...')
    half = GC_WINDOW // 2
    gc_vals = np.full(len(sites_df), np.nan, dtype=np.float32)

    for chrom in sites_df['chrom'].unique():
        try:
            chrom_len = ref.get_reference_length(chrom)
        except KeyError:
            continue
        mask = sites_df['chrom'] == chrom
        positions = sites_df.loc[mask, 'start'].values
        idx = np.where(mask.values)[0]

        for i, pos in zip(idx, positions):
            s = max(0, pos - half)
            e = min(chrom_len, pos + half)
            seq = ref.fetch(chrom, s, e).upper()
            gc = (seq.count('G') + seq.count('C')) / len(seq) if seq else np.nan
            gc_vals[i] = gc

    return gc_vals


def analyze_sample(sample_name, frozen_id, ensilicated_id, islands_df, ref):
    """run full discordance analysis for one sample."""
    print(f'\n=== {sample_name} ===')

    # load and merge
    df_frozen = load_bedgraph(frozen_id)
    df_ensilicated = load_bedgraph(ensilicated_id)

    # filter to autosomes + chrX
    df_frozen = df_frozen[df_frozen['chrom'].isin(AUTOSOMES)]
    df_ensilicated = df_ensilicated[df_ensilicated['chrom'].isin(AUTOSOMES)]

    merged = pd.merge(df_frozen, df_ensilicated, on=['chrom', 'start', 'end'],
                      suffixes=('_frozen', '_cached'))

    # coverage filter
    merged = merged[(merged['coverage_frozen'] >= MIN_COV) &
                    (merged['coverage_ensilicated'] >= MIN_COV)].copy()

    print(f'  {len(merged):,} sites with coverage >= {MIN_COV} in both')

    # compute discordance
    merged['delta'] = (merged['beta_ensilicated'] - merged['beta_frozen']).abs()
    merged['discordant'] = merged['delta'] > DISCORD_THRESH

    n_total = len(merged)
    n_discord = merged['discordant'].sum()
    frac_discord = n_discord / n_total

    print(f'  {n_discord:,} discordant sites ({frac_discord*100:.3f}%)')

    # a. counts
    summary = {
        'sample': sample_name,
        'total_sites': n_total,
        'discordant_sites': n_discord,
        'discordant_fraction': round(frac_discord, 6),
    }

    # b. coverage distributions
    conc = merged[~merged['discordant']]
    disc = merged[merged['discordant']]

    for label, subset in [('concordant', conc), ('discordant', disc)]:
        cov = subset[['coverage_frozen', 'coverage_ensilicated']].values.flatten()
        summary[f'{label}_coverage_median'] = float(np.median(cov))
        summary[f'{label}_coverage_q25'] = float(np.percentile(cov, 25))
        summary[f'{label}_coverage_q75'] = float(np.percentile(cov, 75))

    # c. per-chromosome breakdown
    chrom_stats = []
    for chrom in AUTOSOMES:
        cm = merged[merged['chrom'] == chrom]
        if len(cm) == 0:
            continue
        nd = cm['discordant'].sum()
        chrom_stats.append({
            'chrom': chrom,
            'total': len(cm),
            'discordant': nd,
            'frac_discordant': nd / len(cm),
        })
    chrom_df = pd.DataFrame(chrom_stats)

    # flag enriched chromosomes (> 2x the genome-wide rate)
    chrom_df['enriched'] = chrom_df['frac_discordant'] > 2 * frac_discord
    enriched = chrom_df[chrom_df['enriched']]['chrom'].tolist()
    summary['enriched_chromosomes'] = ','.join(enriched) if enriched else 'none'

    # d. CpG island context
    merged['cpg_context'] = annotate_cpg_context(merged, islands_df)

    context_stats = []
    for ctx in ['island', 'shore', 'shelf', 'open_sea']:
        cm = merged[merged['cpg_context'] == ctx]
        if len(cm) == 0:
            continue
        nd = cm['discordant'].sum()
        context_stats.append({
            'context': ctx,
            'total': len(cm),
            'frac_of_all': len(cm) / n_total,
            'discordant': nd,
            'frac_discordant': nd / len(cm) if len(cm) > 0 else 0,
        })
    context_df = pd.DataFrame(context_stats)

    # e. GC content
    merged['gc'] = compute_gc_content(merged, ref)

    gc_conc = merged.loc[~merged['discordant'], 'gc'].dropna()
    gc_disc = merged.loc[merged['discordant'], 'gc'].dropna()
    summary['gc_concordant_mean'] = float(gc_conc.mean())
    summary['gc_discordant_mean'] = float(gc_disc.mean())

    # save tables
    chrom_df.to_csv(f'{OUT_DIR}/{sample_name}_chrom_breakdown.tsv', sep='\t', index=False)
    context_df.to_csv(f'{OUT_DIR}/{sample_name}_cpg_context.tsv', sep='\t', index=False)

    # scatter plot: frozen vs ensilicated beta, discordant in red
    print('  plotting...')
    fig, ax = plt.subplots(figsize=(5, 5))
    # subsample concordant sites for plotting (too many to render)
    conc_sub = conc.sample(min(500_000, len(conc)), random_state=42)
    ax.scatter(conc_sub['beta_frozen'], conc_sub['beta_ensilicated'],
               s=0.1, alpha=0.05, c='#888888', rasterized=True)
    if len(disc) > 0:
        disc_sub = disc.sample(min(100_000, len(disc)), random_state=42)
        ax.scatter(disc_sub['beta_frozen'], disc_sub['beta_ensilicated'],
                   s=0.3, alpha=0.15, c='#D32F2F', rasterized=True)
    ax.plot([0, 1], [0, 1], 'k--', lw=0.5, alpha=0.5)
    ax.set_xlabel('Frozen beta')
    ax.set_ylabel('Ensilicated beta')
    ax.set_title(f'{sample_name}: Frozen vs ensilicated methylation')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect('equal')
    ax.text(0.05, 0.92, f'{n_discord:,} discordant ({frac_discord*100:.2f}%)',
            transform=ax.transAxes, fontsize=9, color='#D32F2F')
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/{sample_name}_methylation_scatter.png', dpi=200)
    plt.savefig(f'{FIG_DIR}/{sample_name}_methylation_scatter.pdf')
    plt.close()

    # box plot: coverage distributions
    fig, ax = plt.subplots(figsize=(4, 4))
    conc_cov = conc[['coverage_frozen', 'coverage_ensilicated']].values.flatten()
    disc_cov = disc[['coverage_frozen', 'coverage_ensilicated']].values.flatten()
    # cap at 99th percentile for visualization
    cap = np.percentile(np.concatenate([conc_cov, disc_cov]), 99)
    bp = ax.boxplot([conc_cov[conc_cov <= cap], disc_cov[disc_cov <= cap]],
                    labels=['concordant', 'discordant'],
                    patch_artist=True, showfliers=False,
                    medianprops={'color': 'black'})
    bp['boxes'][0].set_facecolor('#B0BEC5')
    bp['boxes'][1].set_facecolor('#EF9A9A')
    ax.set_ylabel('Coverage [reads]')
    ax.set_title(f'{sample_name}: Coverage by discordance')
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    # add median labels
    for i, data in enumerate([conc_cov, disc_cov]):
        med = np.median(data)
        ax.text(i + 1, med + 1, f'{med:.0f}', ha='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/{sample_name}_coverage_boxplot.png', dpi=200)
    plt.savefig(f'{FIG_DIR}/{sample_name}_coverage_boxplot.pdf')
    plt.close()

    return summary


def main():
    islands_df = load_cpg_islands()
    ref = pysam.FastaFile(REF_FASTA)

    all_summaries = []
    for sample_name, ids in SAMPLES.items():
        summary = analyze_sample(sample_name, ids['frozen'], ids['ensilicated'],
                                 islands_df, ref)
        all_summaries.append(summary)

    # save combined summary
    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv(f'{OUT_DIR}/discordance_summary.tsv', sep='\t', index=False)
    print(f'\nsummary saved to {OUT_DIR}/discordance_summary.tsv')
    print(summary_df.to_string(index=False))

    ref.close()


if __name__ == '__main__':
    main()
