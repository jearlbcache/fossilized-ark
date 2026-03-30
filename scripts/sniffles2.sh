#!/usr/bin/env bash
# sniffles2.sh — Call structural variants with Sniffles2
#
# runs Sniffles2 on each aligned BAM to detect structural variants,
# producing compressed and indexed VCF output.
#
# usage:
#   bash scripts/sniffles2.sh <bam_dir> <out_dir> [threads]
#
# arguments:
#   bam_dir Directory containing .aligned.sorted.bam files
#   out_dir Output directory for SV VCFs
#   threads Thread count (default: 64)
#
# example:
#   bash scripts/sniffles2.sh \
#   archive/data_pipeline/12_giab_hg38_aligned \
#   archive/data_pipeline/06_sv \
#   64
#
# dependencies: sniffles
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

bam_dir="$1"
out_dir="$2"
threads="${3:-64}"

[[ -d "$bam_dir" ]] || error "BAM directory not found: $bam_dir"

# dependency check

command -v sniffles &>/dev/null || error "sniffles not found in PATH"

# main

mkdir -p "$out_dir"

count=0
for bam_file in "$bam_dir"/*.aligned.sorted.bam; do
    [[ -f "$bam_file" ]] || continue
    sample=$(basename "$bam_file" .aligned.sorted.bam)

    info "Calling SVs for $sample ..."
    sniffles \
        -i "$bam_file" \
        -v "$out_dir/${sample}.sv.vcf.gz" \
        --threads "$threads"
    (( ++count ))
done

info "Done — $count sample(s) processed to $out_dir"
