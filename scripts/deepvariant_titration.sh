#!/usr/bin/env bash
# deepvariant_titration.sh — Run DeepVariant on titration BAMs (SSD-staged)
#
# calls variants on ensilicated/frozen titration mixtures of HG002, running
# one at a time with BAM+reference staged to NVMe SSD for faster I/O.
#
# usage:
#   bash scripts/deepvariant_titration.sh <project_root> <ref> <titration_dir>
#
# arguments:
#   project_root Absolute path to project root
#   ref Reference FASTA, relative to project_root
#   titration_dir Directory with titration BAMs, relative to project_root
#
# example:
#   bash scripts/deepvariant_titration.sh \
#   /data/jb/project/giab_stanford \
#   reference/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
#   results/titration
#
# dependencies: docker (with GPU support)
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

project_root="$1"
ref="$2"
titration_dir="$3"
shards="${4:-32}"

[[ -d "$project_root" ]]                    || error "Project root not found: $project_root"
[[ -f "$project_root/$ref" ]]               || error "Reference not found: $project_root/$ref"
[[ -d "$project_root/$titration_dir" ]]     || error "Titration dir not found: $project_root/$titration_dir"

# ssd staging directory

SSD_DIR="/home/cache-jb/dv_staging"
mkdir -p "$SSD_DIR"

cleanup_ssd() {
    info "Cleaning SSD staging area..."
    rm -rf "$SSD_DIR"
}
trap cleanup_ssd EXIT

# stage reference to SSD once (shared across all runs)
ref_basename="$(basename "$ref")"
ref_dir="$SSD_DIR/reference"
mkdir -p "$ref_dir"

if [[ ! -f "$ref_dir/$ref_basename" ]]; then
    info "Staging reference to SSD..."
    cp "$project_root/$ref" "$ref_dir/$ref_basename"
    cp "$project_root/$ref.fai" "$ref_dir/$ref_basename.fai" 2>/dev/null || true
    # copy any .dict or other index files
    for ext in .dict .gzi; do
        [[ -f "$project_root/$ref$ext" ]] && cp "$project_root/$ref$ext" "$ref_dir/$ref_basename$ext"
    done
    info "Reference staged."
fi

# dependency check

command -v docker &>/dev/null || error "docker not found in PATH"

# deepVariant runner (SSD-staged, serial)

run_dv() {
    local pct="$1"
    local bam_name="hg002_mix_${pct}pct_ensilicated.bam"
    local vcf_name="hg002_mix_${pct}pct_ensilicated.deepvariant.vcf.gz"
    local src_bam="$project_root/$titration_dir/$bam_name"
    local dst_vcf="$project_root/$titration_dir/$vcf_name"

    if [[ -f "$dst_vcf" ]]; then
        info "SKIP ${pct}% ensilicated — VCF exists"
        return 0
    fi

    [[ -f "$src_bam" ]] || { warn "BAM not found: $src_bam — skipping"; return 1; }

    # clean leftover tmp files from prior failed runs
    rm -f "$src_bam".tmp.*.bam

    # stage BAM + index to SSD
    local ssd_bam="$SSD_DIR/$bam_name"
    local ssd_vcf="$SSD_DIR/$vcf_name"
    info "Staging ${pct}% BAM to SSD..."
    cp "$src_bam" "$ssd_bam"
    cp "$src_bam.bai" "$ssd_bam.bai" 2>/dev/null || \
        cp "${src_bam%.bam}.bai" "${ssd_bam%.bam}.bai" 2>/dev/null || true

    info "START ${pct}% ensilicated (shards=$shards, SSD-staged)"
    sg docker -c "docker run --gpus all \
        -v ${SSD_DIR}:/ssd \
        -v ${ref_dir}:/ref \
        google/deepvariant:1.6.1 \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=ONT_R104 \
        --ref=/ref/${ref_basename} \
        --reads=/ssd/${bam_name} \
        --output_vcf=/ssd/${vcf_name} \
        --num_shards=${shards}"

    # copy VCF back to HDD
    info "Copying VCF back to $dst_vcf..."
    cp "$ssd_vcf" "$dst_vcf"
    cp "$ssd_vcf.tbi" "$dst_vcf.tbi" 2>/dev/null || true
    # copy visual report if generated
    local report="${vcf_name%.vcf.gz}.visual_report.html"
    [[ -f "$SSD_DIR/$report" ]] && cp "$SSD_DIR/$report" "$project_root/$titration_dir/$report"

    # free SSD space for next BAM
    rm -f "$ssd_bam" "$ssd_bam.bai" "${ssd_bam%.bam}.bai" "$ssd_vcf" "$ssd_vcf.tbi" "$SSD_DIR/$report"

    info "DONE  ${pct}% ensilicated"
}

# run serially

pcts=(0 10 25 50 75 90 100)
total=${#pcts[@]}
info "Running DeepVariant on $total titration BAMs (serial, $shards shards, SSD-staged)"

failed=0
for pct in "${pcts[@]}"; do
    run_dv "$pct" || (( ++failed ))
done

if [[ $failed -gt 0 ]]; then
    warn "$failed job(s) failed out of $total"
fi

info "All titration DeepVariant runs complete."
