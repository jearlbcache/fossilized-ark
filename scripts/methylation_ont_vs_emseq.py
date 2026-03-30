#!/usr/bin/env python3
"""compare ONT methylation (frozen and ensilicated) against EMSeq WGBS truth.

for each sample, joins frozen-ONT, ensilicated-ONT, and EMSeq on CpG position,
computes discordance rates and correlations relative to EMSeq, and generates
scatter plots and a summary bar chart.

usage:
    python scripts/methylation_ont_vs_emseq.py

dependencies: pandas, numpy, scipy, matplotlib
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
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
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

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


def load_ont(sample_id):
    """load ONT combined bedgraph: chrom, start, end, beta, coverage."""
    path = f'{MODKIT_DIR}/{sample_id}.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph'
    print(f'  loading {os.path.basename(path)}...')
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'beta', 'coverage'],
                     dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                            'beta': np.float64, 'coverage': np.int32})
    df = df[df['chrom'].isin(AUTOSOMES)]
    return df


def load_emseq(sample_id):
    """load EMSeq converted bedgraph: chrom, start, end, meth_pct, m_reads, u_reads."""
    path = f'{MODKIT_DIR}/{sample_id}.converted.bedGraph'
    print(f'  loading {os.path.basename(path)}...')
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'meth_pct', 'm_reads', 'u_reads'],
                     dtype={'chrom': str, 'start': np.int32, 'end': np.int32,
                            'meth_pct': np.float64, 'm_reads': np.int32, 'u_reads': np.int32})
    df['coverage'] = df['m_reads'] + df['u_reads']
    df['beta'] = df['meth_pct'] / 100.0
    df = df[df['chrom'].isin(AUTOSOMES)]
    return df[['chrom', 'start', 'end', 'beta', 'coverage']]


def analyze_sample(sample_name, ids):
    """run three-way comparison for one sample."""
    print(f'\n=== {sample_name} ===')

    df_frozen = load_ont(ids['frozen'])
    df_ensilicated = load_ont(ids['ensilicated'])
    df_emseq = load_emseq(ids['emseq'])

    # merge all three on position
    m = pd.merge(df_frozen, df_ensilicated, on=['chrom', 'start', 'end'],
                 suffixes=('_frozen', '_cached'))
    m = pd.merge(m, df_emseq[['chrom', 'start', 'end', 'beta', 'coverage']],
                 on=['chrom', 'start', 'end'])
    m.rename(columns={'beta': 'beta_emseq', 'coverage': 'coverage_emseq'}, inplace=True)

    # coverage filters
    m = m[(m['coverage_frozen'] >= MIN_COV_ONT) &
          (m['coverage_ensilicated'] >= MIN_COV_ONT) &
          (m['coverage_emseq'] >= MIN_COV_EMSEQ)].copy()

    n = len(m)
    print(f'  {n:,} sites passing all filters')

    # compute deltas vs EMSeq
    m['delta_frozen'] = np.abs(m['beta_frozen'] - m['beta_emseq'])
    m['delta_ensilicated'] = np.abs(m['beta_ensilicated'] - m['beta_emseq'])

    # discordance rates
    disc_frozen = (m['delta_frozen'] > DISCORD_THRESH).sum()
    disc_ensilicated = (m['delta_ensilicated'] > DISCORD_THRESH).sum()
    rate_frozen = disc_frozen / n
    rate_ensilicated = disc_ensilicated / n

    # median absolute difference
    mad_frozen = m['delta_frozen'].median()
    mad_ensilicated = m['delta_ensilicated'].median()

    # pearson correlations
    r_frozen, _ = pearsonr(m['beta_frozen'], m['beta_emseq'])
    r_ensilicated, _ = pearsonr(m['beta_ensilicated'], m['beta_emseq'])

    print(f'  frozen vs EMSeq: r={r_frozen:.4f}, discord={rate_frozen*100:.2f}%, MAD={mad_frozen:.4f}')
    print(f'  ensilicated vs EMSeq: r={r_ensilicated:.4f}, discord={rate_ensilicated*100:.2f}%, MAD={mad_ensilicated:.4f}')

    result = {
        'sample': sample_name,
        'total_sites': n,
        'r_frozen_emseq': round(r_frozen, 5),
        'r_ensilicated_emseq': round(r_ensilicated, 5),
        'discord_frozen_emseq': int(disc_frozen),
        'discord_rate_frozen': round(rate_frozen, 6),
        'discord_cached_emseq': int(disc_ensilicated),
        'discord_rate_ensilicated': round(rate_ensilicated, 6),
        'mad_frozen_emseq': round(mad_frozen, 5),
        'mad_ensilicated_emseq': round(mad_ensilicated, 5),
    }

    # scatter plots
    print('  plotting...')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # subsample for plotting
    rng = np.random.default_rng(42)
    sub = rng.choice(n, size=min(500_000, n), replace=False)

    for ax, ont_col, label, r_val, disc_rate, mad_val in [
        (ax1, 'beta_frozen', 'frozen', r_frozen, rate_frozen, mad_frozen),
        (ax2, 'beta_ensilicated', 'ensilicated', r_ensilicated, rate_ensilicated, mad_ensilicated),
    ]:
        ax.hexbin(m[ont_col].values[sub], m['beta_emseq'].values[sub],
                  gridsize=150, cmap='Greys', mincnt=1, linewidths=0)
        ax.plot([0, 1], [0, 1], 'r--', lw=0.5, alpha=0.5)
        display_label = 'Ensilicated' if label == 'ensilicated' else 'Frozen'
        ax.set_xlabel(f'{display_label} ONT beta')
        ax.set_ylabel('EMSeq beta')
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect('equal')
        ax.set_title(f'{sample_name}: {display_label} ONT vs EMSeq')
        annot = f'r = {r_val:.4f}\ndiscord = {disc_rate*100:.2f}%\nMAD = {mad_val:.4f}'
        ax.text(0.05, 0.95, annot, transform=ax.transAxes, fontsize=9,
                va='top', fontfamily='monospace',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/{sample_name}_ont_vs_emseq.png', dpi=200)
    plt.savefig(f'{FIG_DIR}/{sample_name}_ont_vs_emseq.pdf')
    plt.close()

    return result


def main():
    all_results = []
    for sample_name, ids in SAMPLES.items():
        r = analyze_sample(sample_name, ids)
        all_results.append(r)

    # summary table
    summary = pd.DataFrame(all_results)
    summary.to_csv(f'{OUT_DIR}/ont_vs_emseq_summary.tsv', sep='\t', index=False)
    print(f'\n{summary.to_string(index=False)}')

    # grouped bar chart: discordance rates frozen vs ensilicated, per sample
    fig, ax = plt.subplots(figsize=(5, 4))
    x = np.arange(len(all_results))
    w = 0.35
    rates_frozen = [r['discord_rate_frozen'] * 100 for r in all_results]
    rates_ensilicated = [r['discord_rate_ensilicated'] * 100 for r in all_results]
    names = [r['sample'] for r in all_results]

    bars1 = ax.bar(x - w / 2, rates_frozen, w, label='Frozen vs EMSeq',
                   color='#555555', alpha=0.85)
    bars2 = ax.bar(x + w / 2, rates_ensilicated, w, label='Ensilicated vs EMSeq',
                   color='#00D8A4', alpha=0.85)

    for bars in [bars1, bars2]:
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
                    f'{bar.get_height():.1f}%', ha='center', fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel('Discordance rate vs EMSeq [%]')
    ax.set_title('ONT vs EMSeq: Frozen and ensilicated', fontweight='bold')
    ax.legend(frameon=False)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.set_ylim(0, max(rates_frozen + rates_ensilicated) * 1.2)

    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/ont_vs_emseq_discordance_barplot.png', dpi=200)
    plt.savefig(f'{FIG_DIR}/ont_vs_emseq_discordance_barplot.pdf')
    plt.close()
    print(f'\nsaved to {FIG_DIR}/ont_vs_emseq_discordance_barplot.png')


if __name__ == '__main__':
    main()
