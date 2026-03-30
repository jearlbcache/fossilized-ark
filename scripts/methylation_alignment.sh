#!/usr/bin/env bash
# methylation_alignment.sh — Align methylation-tagged uBAMs preserving MM/ML
#
# converts Dorado methylation-tagged unaligned BAMs to FASTQ (preserving
# auxiliary tags via -T '*'), aligns with minimap2, and produces sorted,
# indexed BAMs with methylation tags intact.
#
# usage:
#   bash scripts/methylation_alignment.sh <ref> <ubam_dir> <bam_dir> [threads]
#
# arguments:
#   ref Reference FASTA
#   ubam_dir Directory containing unaligned .bam files from Dorado
#   bam_dir Output directory for aligned BAMs
#   threads Thread count (default: 32)
#
# example:
#   bash scripts/methylation_alignment.sh \
#   reference/hg38.fa \
#   archive/data_pipeline/02_methylation_ubam \
#   archive/data_pipeline/10_methylation_aligned \
#   32
#
# dependencies: samtools, minimap2
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
ubam_dir="$2"
bam_dir="$3"
threads="${4:-32}"

[[ -f "$ref" ]]      || error "Reference not found: $ref"
[[ -d "$ubam_dir" ]] || error "uBAM directory not found: $ubam_dir"

# dependency check

for cmd in samtools minimap2; do
    command -v "$cmd" &>/dev/null || error "$cmd not found in PATH"
done

# main

mkdir -p "$bam_dir"

count=0
for ubam in "$ubam_dir"/*.bam; do
    [[ -f "$ubam" ]] || continue
    sample=$(basename "$ubam" .bam)
    out_bam="$bam_dir/${sample}.methyl.aligned.sorted.bam"

    info "Aligning $sample ..."
    samtools fastq -@ "$threads" -T '*' "$ubam" \
        | minimap2 -ax map-ont -t "$threads" -y -K5g "$ref" - \
        | samtools sort -@ "$threads" -O BAM --write-index \
            -o "$out_bam" -
    (( ++count ))
done

info "Done — $count sample(s) aligned to $bam_dir"
