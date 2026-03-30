#!/usr/bin/env bash
# deepvariant.sh — Call small variants with DeepVariant (ONT model, GPU)
#
# runs Google DeepVariant v1.6.1 via Docker on each aligned BAM file,
# using the ONT_R104 model with GPU acceleration.
#
# usage:
#   bash scripts/deepvariant.sh <project_root> <ref> <bam_dir> <out_dir> [shards]
#
# arguments:
#   project_root Absolute path to the project root (mounted into Docker)
#   ref Reference FASTA, relative to project_root
#   bam_dir BAM directory, relative to project_root
#   out_dir Output directory, relative to project_root
#   shards Number of DeepVariant shards (default: 32)
#
# example:
#   bash scripts/deepvariant.sh \
#   /data/jb/project/giab_stanford \
#   reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
#   archive/data_pipeline/12_giab_hg38_aligned \
#   archive/data_pipeline/05_variants \
#   32
#
# dependencies: docker (with GPU support), sudo
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

[[ $# -ge 4 ]] || usage

project_root="$1"
ref="$2"
bam_dir="$3"
out_dir="$4"
shards="${5:-32}"

[[ -d "$project_root" ]]            || error "Project root not found: $project_root"
[[ -f "$project_root/$ref" ]]       || error "Reference not found: $project_root/$ref"
[[ -d "$project_root/$bam_dir" ]]   || error "BAM directory not found: $project_root/$bam_dir"

# dependency check

command -v docker &>/dev/null || error "docker not found in PATH"

# main

mkdir -p "$project_root/$out_dir"

count=0
for bam_file in "$project_root/$bam_dir"/*.nist.aligned.sorted.bam; do
    [[ -f "$bam_file" ]] || continue
    sample=$(basename "$bam_file" .nist.aligned.sorted.bam)

    info "Calling variants for $sample ..."
    sudo docker run --gpus 1 \
        -v "${project_root}:/input" \
        google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=ONT_R104 \
        --ref="/input/${ref}" \
        --reads="/input/${bam_dir}/${sample}.nist.aligned.sorted.bam" \
        --output_vcf="/input/${out_dir}/${sample}.deepvariant.vcf.gz" \
        --output_gvcf="/input/${out_dir}/${sample}.deepvariant.g.vcf.gz" \
        --num_shards="$shards"
    (( ++count ))
done

info "Done — $count sample(s) processed to $out_dir"
