#!/usr/bin/env bash
# fastq_alignment.sh — Align ONT FASTQ reads to a reference with minimap2
#
# aligns each .fastq.gz in the input directory to the given reference genome,
# producing coordinate-sorted and indexed BAM files.
#
# usage:
#   bash scripts/fastq_alignment.sh <ref> <fastq_dir> <bam_dir> [threads]
#
# arguments:
#   ref Reference FASTA (must be indexed)
#   fastq_dir Directory containing .fastq.gz files
#   bam_dir Output directory for aligned BAMs
#   threads Thread count for minimap2 and samtools (default: 64)
#
# example:
#   bash scripts/fastq_alignment.sh \
#   reference/hg38.fa \
#   archive/data_pipeline/03_concatenated_reads \
#   archive/data_pipeline/04_alignment \
#   64
#
# dependencies: minimap2, samtools
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
fastq_dir="$2"
bam_dir="$3"
threads="${4:-64}"

[[ -f "$ref" ]]      || error "Reference not found: $ref"
[[ -d "$fastq_dir" ]] || error "FASTQ directory not found: $fastq_dir"

# dependency check

for cmd in minimap2 samtools; do
    command -v "$cmd" &>/dev/null || error "$cmd not found in PATH"
done

# main

mkdir -p "$bam_dir"

count=0
for fastq in "$fastq_dir"/*.fastq.gz; do
    [[ -f "$fastq" ]] || continue
    sample=$(basename "$fastq" .fastq.gz)
    out_bam="$bam_dir/${sample}.aligned.sorted.bam"

    info "Aligning $sample ..."
    minimap2 -ax map-ont -t "$threads" "$ref" "$fastq" \
        | samtools sort -@ "$threads" -o "$out_bam"
    samtools index -@ "$threads" "$out_bam"
    (( ++count ))
done

info "Done — $count sample(s) aligned to $bam_dir"
