#!/usr/bin/env python3
"""Cross-preservation Mendelian concordance analysis for an Ashkenazi trio.

Tests whether swapping one ensilicated sample into an otherwise frozen trio
(or vice versa) affects Mendelian inheritance consistency. Runs all 5
combinations: all-frozen, all-ensilicated, and 3 single-swap configurations.

Usage:
    python scripts/mendelian_concordance.py <vcf_dir> <output_tsv> <output_figure>

Example:
    python scripts/mendelian_concordance.py \
        archive/data_pipeline/05_variants \
        results/mendelian_cross_preservation.tsv \
        figures/mendelian_cross_preservation.pdf

Dependencies: matplotlib, numpy
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

if len(sys.argv) < 4:
    sys.exit(f"Usage: {sys.argv[0]} <vcf_dir> <output_tsv> <output_figure>")

vcf_dir = sys.argv[1]
output_tsv = sys.argv[2]
output_fig = sys.argv[3]

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 10,
    'pdf.fonttype': 42,
})


FROZEN = {'HG002': 'CD_3033_GIAB', 'HG003': 'CD_3029_GIAB', 'HG004': 'CD_3031_GIAB'}
CACHED = {'HG002': 'CD_3032_Cache', 'HG003': 'CD_3028_Cache', 'HG004': 'CD_3030_Cache'}

# father=HG003, mother=HG004, son=HG002
# all 8 combinations of frozen/ensilicated for the trio
TRIOS = [
    ('FFF', ('HG003', 'Frozen'), ('HG004', 'Frozen'), ('HG002', 'Frozen')),
    ('CFF', ('HG003', 'Ensilicated'), ('HG004', 'Frozen'), ('HG002', 'Frozen')),
    ('FCF', ('HG003', 'Frozen'), ('HG004', 'Ensilicated'), ('HG002', 'Frozen')),
    ('CCF', ('HG003', 'Ensilicated'), ('HG004', 'Ensilicated'), ('HG002', 'Frozen')),
    ('FFC', ('HG003', 'Frozen'), ('HG004', 'Frozen'), ('HG002', 'Ensilicated')),
    ('CFC', ('HG003', 'Ensilicated'), ('HG004', 'Frozen'), ('HG002', 'Ensilicated')),
    ('FCC', ('HG003', 'Frozen'), ('HG004', 'Ensilicated'), ('HG002', 'Ensilicated')),
    ('CCC', ('HG003', 'Ensilicated'), ('HG004', 'Ensilicated'), ('HG002', 'Ensilicated')),
]



def load_snps(vcf_id):
    """Load PASS SNP genotypes from a DeepVariant VCF."""
    path = os.path.join(vcf_dir, f"{vcf_id}.nist.aligned.sorted.full.vcf")
    snps = {}
    with open(path) as f:
        for line in f:
            if line[0] == '#':
                continue
            parts = line.split('\t', 10)
            if len(parts[3]) != 1 or len(parts[4]) != 1:
                continue
            if parts[6] not in ('PASS', '.'):
                continue
            gt_field = parts[9].split(':')[0]
            sep = '/' if '/' in gt_field else '|' if '|' in gt_field else None
            if not sep:
                continue
            try:
                gt = tuple(sorted(int(a) for a in gt_field.split(sep)))
            except ValueError:
                continue
            snps[(parts[0], int(parts[1]))] = gt
    return snps


def mendelian_check(father_snps, mother_snps, child_snps):
    """Compute Mendelian violation rate on sites genotyped in all three."""
    common = set(father_snps) & set(mother_snps) & set(child_snps)
    violations = 0
    for site in common:
        c = child_snps[site]
        f_alleles = set(father_snps[site])
        m_alleles = set(mother_snps[site])
        valid = any(
            tuple(sorted([fa, ma])) == c
            for fa in f_alleles for ma in m_alleles
        )
        if not valid:
            violations += 1
    return len(common), violations



print("Loading VCFs...")
vcfs = {}
for sample in ['HG002', 'HG003', 'HG004']:
    for cond, mapping in [('Frozen', FROZEN), ('Ensilicated', CACHED)]:
        print(f"  {sample} {cond}...")
        vcfs[(sample, cond)] = load_snps(mapping[sample])


print()
results = []
for name, father_key, mother_key, child_key in TRIOS:
    n_sites, n_viol = mendelian_check(vcfs[father_key], vcfs[mother_key], vcfs[child_key])
    rate = n_viol / n_sites * 100
    results.append({
        'trio': name,
        'father_sample': father_key[0], 'father_condition': father_key[1],
        'mother_sample': mother_key[0], 'mother_condition': mother_key[1],
        'son_sample': child_key[0], 'son_condition': child_key[1],
        'sites_genotyped': n_sites, 'violations': n_viol,
        'violation_rate_pct': round(rate, 4),
    })
    print(f"  {name:<15s}: {n_viol:,} / {n_sites:,} = {rate:.4f}%")


os.makedirs(os.path.dirname(output_tsv), exist_ok=True)
header = list(results[0].keys())
with open(output_tsv, 'w') as f:
    f.write('\t'.join(header) + '\n')
    for row in results:
        f.write('\t'.join(str(row[k]) for k in header) + '\n')
print(f"\nSaved: {output_tsv}")


os.makedirs(os.path.dirname(output_fig), exist_ok=True)
import seaborn as sns

# build lookup: (father_cond, mother_cond, son_cond) -> violation rate
lookup = {}
for r in results:
    key = (r['father_condition'], r['mother_condition'], r['son_condition'])
    lookup[key] = r['violation_rate_pct']

# expand to all 8 combinations by running any missing ones
parents = ['Frozen', 'Ensilicated']
labels = ['Frozen', 'Ensilicated']

fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharey=True)

for ax, son_cond in zip(axes, ['Frozen', 'Ensilicated']):
    matrix = np.full((2, 2), np.nan)
    annot = np.empty((2, 2), dtype=object)
    for i, father_cond in enumerate(parents):
        for j, mother_cond in enumerate(parents):
            key = (father_cond, mother_cond, son_cond)
            val = lookup.get(key, np.nan)
            matrix[i, j] = val
            annot[i, j] = f'{val:.3f}%' if not np.isnan(val) else ''

    sns.heatmap(matrix, annot=annot, fmt='', cmap='RdYlGn_r',
                vmin=0.15, vmax=0.20,
                xticklabels=labels, yticklabels=labels,
                linewidths=1.5, linecolor='white',
                cbar=ax == axes[1],
                cbar_kws={'label': 'Violation rate [%]'} if ax == axes[1] else {},
                ax=ax)
    son_label = 'ensilicated' if son_cond == 'Ensilicated' else 'frozen'
    ax.set_title(f'Son: {son_label}', fontweight='bold')
    ax.set_xlabel('Mother')
    if ax == axes[0]:
        ax.set_ylabel('Father')
    else:
        ax.set_ylabel('')

fig.suptitle('Mendelian concordance across preservation methods',
             fontweight='bold', fontsize=11, y=1.02)
plt.tight_layout()

fig_base = output_fig.rsplit('.', 1)[0]
plt.savefig(f'{fig_base}.pdf', bbox_inches='tight')
plt.savefig(f'{fig_base}.png', dpi=300, bbox_inches='tight')
print(f"Saved: {fig_base}.pdf/.png")
