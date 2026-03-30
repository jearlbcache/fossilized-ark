#!/usr/bin/env bash
# merge_ubam.sh — Merge per-read unaligned BAMs from Dorado into one per sample
#
# finds bam_pass/ directories under each sample folder and merges all BAM
# files into a single unaligned BAM using samtools cat.
#
# usage:
#   bash scripts/merge_ubam.sh <source_dir> <dest_dir> [sample_glob] [threads]
#
# arguments:
#   source_dir Root directory containing per-sample Dorado output folders
#   dest_dir Output directory for merged BAMs
#   sample_glob Glob pattern for sample folders (default: CD_*)
#   threads Samtools thread count (default: 32)
#
# example:
#   bash scripts/merge_ubam.sh \
#   /media/cache-jb/ExternalHDD/Cache_Fossilized_Ark_Batch1 \
#   archive/data_pipeline/02_methylation_ubam \
#   "CD_30*" 32
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

[[ $# -ge 2 ]] || usage

source_dir="$1"
dest_dir="$2"
sample_glob="${3:-CD_*}"
threads="${4:-32}"

[[ -d "$source_dir" ]] || error "Source directory not found: $source_dir"

# dependency check

command -v samtools &>/dev/null || error "samtools not found in PATH"

# main

mkdir -p "$dest_dir"

count=0
for parent_dir in "$source_dir"/$sample_glob; do
    [[ -d "$parent_dir" ]] || continue

    sample=$(basename "$parent_dir")
    bam_pass_dir=$(find "$parent_dir" -type d -name "bam_pass" -print -quit)

    if [[ -z "$bam_pass_dir" || ! -d "$bam_pass_dir" ]]; then
        warn "No bam_pass/ directory found under $parent_dir — skipping"
        continue
    fi

    merged_bam="${dest_dir}/${sample}.bam"
    info "Merging $sample ..."
    samtools cat -@ "$threads" -o "$merged_bam" "$bam_pass_dir"/*.bam
    samtools index -@ "$threads" "$merged_bam"
    (( ++count ))
done

info "Done — $count sample(s) merged to $dest_dir"
