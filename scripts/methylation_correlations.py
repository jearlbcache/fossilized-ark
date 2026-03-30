#!/usr/bin/env python3
"""Compute pairwise Pearson correlations between methylation bedGraph files.

Usage: python methylation_correlations.py <bed_dir> <output_csv>
Example: python methylation_correlations.py archive/data_pipeline/11_modkit results/all_pairs_correlations.csv
"""
import os
import sys
import itertools
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm

bed_dir = sys.argv[1] if len(sys.argv) > 1 else sys.exit(f"Usage: {sys.argv[0]} <bed_dir> <output_csv>")
output_csv = sys.argv[2]


def load_bedgraph(filepath):
    """Load a bedGraph file, returning a DataFrame with chrom/start/end/meth_percent."""
    fname = os.path.basename(filepath)
    if "converted.bedGraph" in fname:
        df = pd.read_csv(filepath, sep="\t", header=None,
                         names=["chrom", "start", "end", "meth_percent", "m_reads", "u_reads"],
                         usecols=["chrom", "start", "end", "meth_percent"])
    elif "m_CG0_combined" in fname:
        df = pd.read_csv(filepath, sep="\t", header=None,
                         names=["chrom", "start", "end", "meth_frac", "coverage"],
                         usecols=["chrom", "start", "end", "meth_frac"])
        df["meth_percent"] = df.pop("meth_frac") * 100
    else:
        return None
    return df


def main():
    all_files = [
        os.path.join(bed_dir, f)
        for f in os.listdir(bed_dir)
        if f.endswith((".bedgraph", ".bedGraph"))
        and "h_CG0_combined" not in f
        and ("converted.bedGraph" in f or "m_CG0_combined" in f)
    ]

    # cache parsed files to avoid re-reading per pair
    cache = {}
    for f in all_files:
        cache[f] = load_bedgraph(f)

    # use combinations (not permutations) — correlation is symmetric
    results = []
    for fA, fB in tqdm(list(itertools.combinations(all_files, 2)),
                       desc="Computing correlations", unit="pair"):
        dfA, dfB = cache[fA], cache[fB]
        if dfA is None or dfB is None:
            continue

        merged = pd.merge(dfA, dfB, on=["chrom", "start", "end"], suffixes=(".A", ".B"))
        if merged.empty:
            r_val, p_val, overlap = float("nan"), float("nan"), 0
        else:
            r_val, p_val = pearsonr(merged["meth_percent.A"], merged["meth_percent.B"])
            overlap = len(merged)

        results.append({
            "FileA": os.path.basename(fA),
            "FileB": os.path.basename(fB),
            "pearson_r": r_val,
            "p_value": p_val,
            "n_overlap": overlap,
        })

    res_df = pd.DataFrame(results)
    res_df.to_csv(output_csv, index=False)
    print(f"Saved {len(res_df)} pairs to {output_csv}")


if __name__ == "__main__":
    main()
