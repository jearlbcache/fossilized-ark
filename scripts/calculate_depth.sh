#!/usr/bin/env bash
# calculate_depth.sh — Calculate per-chromosome mean coverage
#
# uses samtools coverage to compute mean depth per chromosome for each BAM
# in the given directory. Outputs a TSV with columns: chromosome, mean_depth.
#
# usage:
#   bash scripts/calculate_depth.sh <bam_dir> [out_dir]
#
# arguments:
#   bam_dir Directory containing .bam files
#   out_dir Output directory for coverage TSVs (default: same as bam_dir)
#
# example:
#   bash scripts/calculate_depth.sh \
#   archive/data_pipeline/12_giab_hg38_aligned \
#   results/coverage
#
# dependencies: samtools
set -euo pipefail

# logging

readonly SCRIPT_NAME="$(basename "$0")"

info()  { printf '[%s] %s\n' "$SCRIPT_NAME" "$*"; }
warn()  { printf '[%s] WARN: %s\n' "$SCRIPT_NAME" "$*" >&2; }
error() { printf '[%s] ERROR: %s\n' "$SCRIPT_NAME" "$*" >&2; exit 1; }

# arguments

usage() {
    sed -n '/^# Usage:/,/^# Dependencies:/{ /^# Dependencies:/d; s/^# //; s/^#//; p }' "$0"
    exit 1
}

[[ $# -ge 1 ]] || usage

bam_dir="$1"
out_dir="${2:-$bam_dir}"

[[ -d "$bam_dir" ]] || error "BAM directory not found: $bam_dir"

# dependency check

command -v samtools &>/dev/null || error "samtools not found in PATH"

# main

mkdir -p "$out_dir"

count=0
for bam_file in "$bam_dir"/*.bam; do
    [[ -f "$bam_file" ]] || continue
    sample=$(basename "$bam_file" .bam)
    output_file="$out_dir/${sample}_chr_coverage.tsv"

    info "Computing coverage for $sample ..."
    printf 'chromosome\tmean_depth\n' > "$output_file"
    samtools coverage "$bam_file" \
        | awk -F'\t' 'NR > 1 { print $1, $7 }' OFS='\t' \
        >> "$output_file"
    (( ++count ))
done

info "Done — $count sample(s) processed to $out_dir"
