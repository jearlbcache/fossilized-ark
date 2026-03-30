#!/usr/bin/env python3
"""investigate EMSeq vs ONT discordance rate differences across samples.

checks whether the higher discordance in HG003/HG004 vs HG002 is explained
by coverage differences, by stratifying discordance by ONT coverage bins.

usage:
    python scripts/methylation_coverage_investigation.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

PROJECT = '/data/jb/project/giab_stanford'
MODKIT_DIR = f'{PROJECT}/archive/data_pipeline/11_modkit'
OUT_DIR = f'{PROJECT}/results/methylation_discordance'
FIG_DIR = f'{PROJECT}/figures'

SAMPLES = {
    'HG002': {
        'frozen': 'CD_3033_GIAB',
        'ensilicated': 'CD_3032_Cache',
        'emseq': 'EMSeq_HG002_LAB01_REP01',
    },
    'HG003': {
        'frozen': 'CD_3031_GIAB',
        'ensilicated': 'CD_3030_Cache',
        'emseq': 'EMSeq_HG003_LAB01_REP01',
    },
    'HG004': {
        'frozen': 'CD_3029_GIAB',
        'ensilicated': 'CD_3028_Cache',
        'emseq': 'EMSeq_HG004_LAB01_REP01',
    },
}

MIN_COV_ONT = 10
MIN_COV_EMSEQ = 10
DISCORD_THRESH = 0.20
AUTOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX']
COV_BINS = [(10, 15), (15, 20), (20, 25), (25, 30), (30, 999)]
COV_BIN_LABELS = ['10-15x', '15-20x', '20-25x', '25-30x', '30x+']


def load_ont(sample_id):
    path = f'{MODKIT_DIR}/{sample_id}.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph'
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'beta', 'coverage'],
                     dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                            'beta': np.float64, 'coverage': np.int32})
    return df[df['chrom'].isin(AUTOSOMES)]


def load_emseq(sample_id):
    path = f'{MODKIT_DIR}/{sample_id}.converted.bedGraph'
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'meth_pct', 'm_reads', 'u_reads'],
                     dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                            'meth_pct': np.float64, 'm_reads': np.int32, 'u_reads': np.int32})
    df['coverage'] = df['m_reads'] + df['u_reads']
    df['beta'] = df['meth_pct'] / 100.0
    return df[df['chrom'].isin(AUTOSOMES)]


def coverage_percentiles(cov, label):
    pcts = [10, 25, 50, 75, 90]
    vals = np.percentile(cov, pcts)
    return {f'{label}_p{p}': v for p, v in zip(pcts, vals)}


def main():
    print('=' * 70)
    print('1. EMSeq file metadata')
    print('=' * 70)
    for name, ids in SAMPLES.items():
        path = f'{MODKIT_DIR}/{ids["emseq"]}.converted.bedGraph'
        print(f'\n  {name}: {os.path.basename(path)}')
        # all EMSeq files have same LAB01_REP01 suffix
        print(f'    lab: LAB01, replicate: REP01 (from filename)')

    print('\n' + '=' * 70)
    print('2. coverage stats per dataset')
    print('=' * 70)

    all_stats = []

    for name, ids in SAMPLES.items():
        print(f'\n--- {name} ---')
        stats = {'sample': name}

        # EMSeq
        df_e = load_emseq(ids['emseq'])
        stats['emseq_total_sites'] = len(df_e)
        stats['emseq_sites_cov10'] = (df_e['coverage'] >= MIN_COV_EMSEQ).sum()
        pcts = coverage_percentiles(df_e['coverage'].values, 'emseq_cov')
        stats.update(pcts)
        print(f'  EMSeq: {len(df_e):,} sites, median cov = {pcts["emseq_cov_p50"]:.0f}')

        # frozen ONT
        df_f = load_ont(ids['frozen'])
        stats['ont_frozen_total_sites'] = len(df_f)
        pcts_f = coverage_percentiles(df_f['coverage'].values, 'ont_frozen_cov')
        stats.update(pcts_f)
        print(f'  ONT frozen: {len(df_f):,} sites, median cov = {pcts_f["ont_frozen_cov_p50"]:.0f}')

        # ensilicated ONT
        df_c = load_ont(ids['ensilicated'])
        stats['ont_ensilicated_total_sites'] = len(df_c)
        pcts_c = coverage_percentiles(df_c['coverage'].values, 'ont_ensilicated_cov')
        stats.update(pcts_c)
        print(f'  ONT ensilicated: {len(df_c):,} sites, median cov = {pcts_c["ont_ensilicated_cov_p50"]:.0f}')

        # merge all three
        m = pd.merge(df_f, df_c, on=['chrom', 'start', 'end'], suffixes=('_frozen', '_cached'))
        m = pd.merge(m, df_e[['chrom', 'start', 'end', 'beta', 'coverage']],
                     on=['chrom', 'start', 'end'])
        m.rename(columns={'beta': 'beta_emseq', 'coverage': 'coverage_emseq'}, inplace=True)
        m = m[(m['coverage_frozen'] >= MIN_COV_ONT) &
              (m['coverage_ensilicated'] >= MIN_COV_ONT) &
              (m['coverage_emseq'] >= MIN_COV_EMSEQ)].copy()

        stats['sites_passing_filters'] = len(m)
        print(f'  sites passing all filters: {len(m):,}')

        # discordance stratified by ONT coverage
        m['min_ont_cov'] = np.minimum(m['coverage_frozen'], m['coverage_ensilicated'])
        m['delta_frozen'] = np.abs(m['beta_frozen'] - m['beta_emseq'])
        m['delta_ensilicated'] = np.abs(m['beta_ensilicated'] - m['beta_emseq'])

        print(f'\n  discordance stratified by min ONT coverage:')
        print(f'  {"bin":<10s}  {"n_sites":>10s}  {"frozen_disc":>12s}  {"ensilicated_disc":>12s}')
        for (lo, hi), label in zip(COV_BINS, COV_BIN_LABELS):
            mask = (m['min_ont_cov'] >= lo) & (m['min_ont_cov'] < hi)
            ns = mask.sum()
            if ns == 0:
                continue
            df_frozen = (m.loc[mask, 'delta_frozen'] > DISCORD_THRESH).mean()
            dc_ensilicated = (m.loc[mask, 'delta_ensilicated'] > DISCORD_THRESH).mean()
            print(f'  {label:<10s}  {ns:>10,}  {df_frozen*100:>11.2f}%  {dc_ensilicated*100:>11.2f}%')
            stats[f'n_{label}'] = int(ns)
            stats[f'frozen_disc_{label}'] = round(df_frozen, 6)
            stats[f'ensilicated_disc_{label}'] = round(dc_ensilicated, 6)

        all_stats.append(stats)

    # save summary
    summary = pd.DataFrame(all_stats)
    summary.to_csv(f'{OUT_DIR}/coverage_investigation.tsv', sep='\t', index=False)

    # stratified discordance figure
    print('\n\nplotting stratified discordance...')
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharey=True)

    for ax, stats in zip(axes, all_stats):
        name = stats['sample']
        x = np.arange(len(COV_BIN_LABELS))
        w = 0.35
        frozen_rates = [stats.get(f'frozen_disc_{l}', 0) * 100 for l in COV_BIN_LABELS]
        cached_rates = [stats.get(f'ensilicated_disc_{l}', 0) * 100 for l in COV_BIN_LABELS]

        ax.bar(x - w / 2, frozen_rates, w, label='Frozen vs EMSeq', color='#555555', alpha=0.85)
        ax.bar(x + w / 2, cached_rates, w, label='Ensilicated vs EMSeq', color='#00D8A4', alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels(COV_BIN_LABELS, rotation=30, ha='right')
        ax.set_title(name, fontweight='bold')
        if ax == axes[0]:
            ax.set_ylabel('Discordance rate vs EMSeq [%]')
            ax.legend(frameon=False, fontsize=8)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    fig.suptitle('Discordance stratified by ONT coverage', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/ont_vs_emseq_stratified_coverage.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/ont_vs_emseq_stratified_coverage.pdf', bbox_inches='tight')
    plt.close()
    print(f'saved to {FIG_DIR}/ont_vs_emseq_stratified_coverage.png')


if __name__ == '__main__':
    main()
