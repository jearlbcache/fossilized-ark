#!/usr/bin/env bash
# gc.sh — Compute GC bias profiles (coverage vs GC content)
#
# creates fixed-size genomic windows, computes GC content per window using
# bedtools nuc, then calculates mean read coverage per window for each BAM.
# outputs a merged table of coverage + GC per sample.
#
# usage:
#   bash scripts/gc.sh <ref> <bam_dir> <out_dir> [window_size]
#
# arguments:
#   ref Reference FASTA (will be indexed if .fai is missing)
#   bam_dir Directory containing aligned .bam files
#   out_dir Output directory for GC bias tables
#   window_size Genomic window size in bp (default: 100)
#
# example:
#   bash scripts/gc.sh \
#   reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
#   archive/data_pipeline/12_giab_hg38_aligned \
#   archive/data_pipeline/20_gc \
#   100
#
# dependencies: samtools, bedtools
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
bam_dir="$2"
out_dir="$3"
window_size="${4:-100}"

[[ -f "$ref" ]]     || error "Reference not found: $ref"
[[ -d "$bam_dir" ]] || error "BAM directory not found: $bam_dir"

# dependency check

for cmd in samtools bedtools; do
    command -v "$cmd" &>/dev/null || error "$cmd not found in PATH"
done

# sample map

declare -A SAMPLES=(
    [HG002_frozen]="$bam_dir/CD_3033_GIAB.nist.aligned.sorted.bam"
    [HG002_ensilicated]="$bam_dir/CD_3032_Cache.nist.aligned.sorted.bam"
    [HG003_frozen]="$bam_dir/CD_3031_GIAB.nist.aligned.sorted.bam"
    [HG003_ensilicated]="$bam_dir/CD_3030_Cache.nist.aligned.sorted.bam"
    [HG004_frozen]="$bam_dir/CD_3029_GIAB.nist.aligned.sorted.bam"
    [HG004_ensilicated]="$bam_dir/CD_3028_Cache.nist.aligned.sorted.bam"
)

# prepare windows

mkdir -p "$out_dir"

if [[ ! -f "${ref}.fai" ]]; then
    info "Indexing reference ..."
    samtools faidx "$ref"
fi

genome_file="$out_dir/genome_file.txt"
cut -f1,2 "${ref}.fai" > "$genome_file"

windows_bed="$out_dir/windows_${window_size}.bed"
gc_bed="$out_dir/windows_${window_size}_gc.bed"

info "Creating ${window_size}bp windows and computing GC content ..."
bedtools makewindows -g "$genome_file" -w "$window_size" > "$windows_bed"
bedtools nuc -fi "$ref" -bed "$windows_bed" \
    | grep -v '^#' \
    | sort -k1,1 -k2,2n \
    > "$gc_bed"

# compute per-sample coverage

count=0
for label in "${!SAMPLES[@]}"; do
    bam="${SAMPLES[$label]}"
    [[ -f "$bam" ]] || { warn "BAM not found: $bam — skipping $label"; continue; }

    info "Computing coverage for $label ..."
    bedtools coverage -a "$windows_bed" -b "$bam" -mean -sorted -g "$genome_file" \
        | sort -k1,1 -k2,2n \
        | paste - "$gc_bed" \
        > "$out_dir/coverage_gc_${label}.txt"
    (( ++count ))
done

info "Done — $count sample(s) processed to $out_dir"
