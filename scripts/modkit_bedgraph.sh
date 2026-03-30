#!/usr/bin/env bash
# modkit_bedgraph.sh — Extract CpG methylation from BAMs with modkit
#
# runs modkit pileup on each aligned BAM to produce CpG methylation
# bedGraph files with combined-strand calls.
#
# usage:
#   bash scripts/modkit_bedgraph.sh <ref> <input_dir> <output_dir> [threads]
#
# arguments:
#   ref Reference FASTA
#   input_dir Directory containing methylation-aligned .bam files
#   output_dir Output directory for bedGraph files
#   threads Thread count (default: 192)
#
# example:
#   bash scripts/modkit_bedgraph.sh \
#   reference/hg38.fa \
#   archive/data_pipeline/10_methylation_aligned \
#   archive/data_pipeline/11_modkit \
#   192
#
# dependencies: modkit
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

[[ $# -ge 3 ]] || usage

ref="$1"
input_dir="$2"
output_dir="$3"
threads="${4:-192}"

[[ -f "$ref" ]]       || error "Reference not found: $ref"
[[ -d "$input_dir" ]] || error "Input directory not found: $input_dir"

# dependency check

command -v modkit &>/dev/null || error "modkit not found in PATH"

# main

mkdir -p "$output_dir"

count=0
for bamfile in "$input_dir"/*.bam; do
    [[ -f "$bamfile" ]] || continue
    prefix=$(basename "$bamfile" .bam)

    info "Processing $prefix ..."
    modkit pileup "$bamfile" "$output_dir" \
        --threads "$threads" \
        --ref "$ref" \
        --combine-strands \
        --cpg \
        --bedgraph \
        --prefix "${prefix}.cpg"
    (( ++count ))
done

info "Done — $count sample(s) processed to $output_dir"
