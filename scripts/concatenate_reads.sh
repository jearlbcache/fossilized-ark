#!/usr/bin/env bash
# concatenate_reads.sh — Concatenate FASTQ files from ONT flow cell runs
#
# finds fastq_pass/ directories under each sample folder and concatenates
# all .fastq.gz files into a single FASTQ per sample.
#
# usage:
#   bash scripts/concatenate_reads.sh <source_dir> <dest_dir> [sample_glob]
#
# arguments:
#   source_dir Root directory containing per-sample run folders
#   dest_dir Output directory for concatenated FASTQs
#   sample_glob Glob pattern for sample folders (default: CD_*)
#
# example:
#   bash scripts/concatenate_reads.sh \
#   /media/cache-jb/ExternalHDD/Cache_Fossilized_Ark_Batch2 \
#   archive/data_pipeline/13_concatenated_reads_batch2 \
#   "CD_38*"
#
# dependencies: none (coreutils only)
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

[[ -d "$source_dir" ]] || error "Source directory not found: $source_dir"

# main

mkdir -p "$dest_dir"

count=0
for parent_dir in "$source_dir"/$sample_glob; do
    [[ -d "$parent_dir" ]] || continue

    sample=$(basename "$parent_dir")
    fastq_pass_dir=$(find "$parent_dir" -type d -name "fastq_pass" -print -quit)

    if [[ -z "$fastq_pass_dir" || ! -d "$fastq_pass_dir" ]]; then
        warn "No fastq_pass/ directory found under $parent_dir — skipping"
        continue
    fi

    info "Concatenating $sample ..."
    cat "$fastq_pass_dir"/*.fastq.gz > "$dest_dir/${sample}.fastq.gz"
    (( ++count ))
done

info "Done — $count sample(s) concatenated to $dest_dir"
