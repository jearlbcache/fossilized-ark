#!/usr/bin/env python3
"""binomial null model for methylation discordance.

tests whether observed discordance between frozen and ensilicated methylation
can be explained by sampling noise alone. for each CpG site, models both
observations as binomial draws from a shared true methylation rate, then
computes the expected number of discordant sites (|delta beta| > 0.20).

usage:
    python scripts/methylation_binomial_null.py

dependencies: pandas, numpy, scipy, matplotlib
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})

# configuration
PROJECT = '/data/jb/project/giab_stanford'
MODKIT_DIR = f'{PROJECT}/archive/data_pipeline/11_modkit'
OUT_DIR = f'{PROJECT}/results/methylation_discordance/binomial_null'
FIG_DIR = f'{PROJECT}/figures'
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = {
    'HG002': {'frozen': 'CD_3033_GIAB', 'ensilicated': 'CD_3032_Cache'},
    'HG003': {'frozen': 'CD_3031_GIAB', 'ensilicated': 'CD_3030_Cache'},
    'HG004': {'frozen': 'CD_3029_GIAB', 'ensilicated': 'CD_3028_Cache'},
}

MIN_COV = 10
DISCORD_THRESH = 0.20
AUTOSOMES = [f'chr{i}' for i in range(1, 23)] + ['chrX']


def load_and_merge(frozen_id, ensilicated_id):
    """load bedgraphs, merge on position, filter by coverage."""
    cols = ['chrom', 'start', 'end', 'beta', 'coverage']
    dtypes = {'chrom': str, 'start': np.int32, 'end': np.int32,
              'beta': np.float64, 'coverage': np.int32}

    path_f = f'{MODKIT_DIR}/{frozen_id}.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph'
    path_c = f'{MODKIT_DIR}/{ensilicated_id}.methyl.aligned.sorted.cpg_m_CG0_combined.bedgraph'

    df_f = pd.read_csv(path_f, sep='\t', header=None, names=cols, dtype=dtypes)
    df_c = pd.read_csv(path_c, sep='\t', header=None, names=cols, dtype=dtypes)

    df_f = df_f[df_f['chrom'].isin(AUTOSOMES)]
    df_c = df_c[df_c['chrom'].isin(AUTOSOMES)]

    m = pd.merge(df_f, df_c, on=['chrom', 'start', 'end'], suffixes=('_frozen', '_cached'))
    m = m[(m['coverage_frozen'] >= MIN_COV) & (m['coverage_ensilicated'] >= MIN_COV)]
    return m


def analyze_sample(sample_name, frozen_id, ensilicated_id):
    """run binomial null model analysis for one sample."""
    print(f'\n=== {sample_name} ===')
    print(f'  loading and merging...')
    df = load_and_merge(frozen_id, ensilicated_id)
    n = len(df)
    print(f'  {n:,} sites')

    beta_f = df['beta_frozen'].values
    beta_c = df['beta_ensilicated'].values
    n_f = df['coverage_frozen'].values.astype(np.float64)
    n_c = df['coverage_ensilicated'].values.astype(np.float64)

    # observed
    delta = np.abs(beta_f - beta_c)
    obs_discord = delta > DISCORD_THRESH
    n_obs = obs_discord.sum()

    # binomial null model
    # coverage-weighted pooled estimate of true methylation
    p_hat = (beta_f * n_f + beta_c * n_c) / (n_f + n_c)

    # expected variance of (beta_f - beta_c) under independent binomial sampling
    var_diff = p_hat * (1 - p_hat) / n_f + p_hat * (1 - p_hat) / n_c

    # avoid division by zero for p_hat = 0 or 1
    sd_diff = np.sqrt(np.maximum(var_diff, 1e-20))

    # P(|delta| > 0.20) per site under normal approximation
    p_discord = 2 * (1 - norm.cdf(DISCORD_THRESH / sd_diff))

    n_expected = p_discord.sum()

    print(f'  observed discordant:  {n_obs:>10,} ({n_obs/n*100:.3f}%)')
    print(f'  expected discordant:  {n_expected:>10,.0f} ({n_expected/n*100:.3f}%)')
    print(f'  ratio obs/exp:        {n_obs/n_expected:.2f}')

    result = {
        'sample': sample_name,
        'total_sites': n,
        'observed_discordant': int(n_obs),
        'observed_rate': round(n_obs / n, 6),
        'expected_discordant': round(n_expected, 0),
        'expected_rate': round(n_expected / n, 6),
        'obs_over_exp': round(n_obs / n_expected, 3),
    }

    # panel A: discordance rate vs coverage (observed and expected)
    print('  plotting panel A (discordance vs coverage)...')
    min_cov = df[['coverage_frozen', 'coverage_ensilicated']].min(axis=1).values
    max_cov_bin = int(np.percentile(min_cov, 99))
    bins = np.arange(MIN_COV, max_cov_bin + 5, 5)

    bin_idx = np.digitize(min_cov, bins) - 1
    obs_rates = []
    exp_rates = []
    bin_centers = []

    for i in range(len(bins) - 1):
        mask = bin_idx == i
        if mask.sum() < 100:
            continue
        bin_centers.append((bins[i] + bins[i + 1]) / 2)
        obs_rates.append(obs_discord[mask].mean())
        exp_rates.append(p_discord[mask].mean())

    # panel B: histogram of observed |delta beta| vs expected distribution
    print('  plotting panel B (delta beta distribution)...')
    # expected distribution: average of folded-normal densities across sites
    # for each site, the delta follows N(0, var_diff), so |delta| is folded-normal
    # aggregate density at a grid of delta values
    delta_grid = np.linspace(0, 0.5, 500)
    # vectorized: for each grid point, average the folded-normal density across sites
    # folded-normal pdf at x for N(0, sigma^2) = 2 * phi(x / sigma) / sigma for x >= 0
    # use a random subsample for speed
    rng = np.random.default_rng(42)
    sub_idx = rng.choice(len(sd_diff), size=min(2_000_000, len(sd_diff)), replace=False)
    sd_sub = sd_diff[sub_idx]

    # compute expected density: average of 2*phi(x/sigma)/sigma across sites
    # do this in chunks to avoid memory explosion
    chunk_size = 200_000
    exp_density = np.zeros(len(delta_grid))
    n_chunks = (len(sd_sub) + chunk_size - 1) // chunk_size
    for ci in range(n_chunks):
        s = ci * chunk_size
        e = min(s + chunk_size, len(sd_sub))
        sds = sd_sub[s:e, np.newaxis]  # [chunk, 1]
        grid = delta_grid[np.newaxis, :]  # [1, grid]
        # folded normal density
        density = 2 * norm.pdf(grid / sds) / sds  # [chunk, grid]
        exp_density += density.sum(axis=0)
    exp_density /= len(sd_sub)

    # plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(bin_centers, obs_rates, 'o-', color='#D32F2F', markersize=4, label='Observed')
    ax1.plot(bin_centers, exp_rates, 's--', color='#1565C0', markersize=4, label='Expected [binomial null]')
    ax1.set_xlabel('Minimum coverage [reads]')
    ax1.set_ylabel('Fraction discordant')
    ax1.set_title(f'{sample_name}: Discordance vs coverage')
    ax1.legend(frameon=False)
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)

    tail_mask = delta >= 0.05
    ax2.hist(delta[tail_mask], bins=225, range=(0.05, 0.5), density=True, alpha=0.5,
             color='#D32F2F', label='Observed')
    grid_mask = delta_grid >= 0.05
    ax2.plot(delta_grid[grid_mask], exp_density[grid_mask], color='#1565C0', lw=2,
             label='Expected [binomial null]')
    ax2.axvline(DISCORD_THRESH, color='gray', ls='--', lw=0.8, alpha=0.7)
    ax2.set_xlim(0.05, 0.50)
    ax2.set_xlabel('|Delta beta|')
    ax2.set_ylabel('Density')
    ax2.set_title(f'{sample_name}: Delta beta distribution')
    ax2.legend(frameon=False)
    obs_exp = n_obs / n_expected
    ax2.text(0.95, 0.75, f'obs/exp = {obs_exp:.2f}',
             transform=ax2.transAxes, ha='right', fontsize=10, fontweight='bold')
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/{sample_name}_binomial_null.png', dpi=200)
    plt.savefig(f'{FIG_DIR}/{sample_name}_binomial_null.pdf')
    plt.close()
    print(f'  saved {FIG_DIR}/{sample_name}_binomial_null.png')

    return result


def main():
    all_results = []
    for sample_name, ids in SAMPLES.items():
        r = analyze_sample(sample_name, ids['frozen'], ids['ensilicated'])
        all_results.append(r)

    summary = pd.DataFrame(all_results)
    summary.to_csv(f'{OUT_DIR}/binomial_null_summary.tsv', sep='\t', index=False)
    print(f'\n{summary.to_string(index=False)}')
    print(f'\nsaved to {OUT_DIR}/binomial_null_summary.tsv')


if __name__ == '__main__':
    main()
