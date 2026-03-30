#!/usr/bin/env bash
# fix_vcf.sh — Resolve symbolic SV alleles and normalize a VCF
#
# replaces <DEL> with the reference sequence, drops unresolved symbolic
# alleles (<INV>, <DUP>, <INS>), and runs bcftools norm for left-alignment
# and multiallelic splitting. Output is bgzipped and tabix-indexed.
#
# usage:
#   bash scripts/fix_vcf.sh <ref> <in_vcf> <out_vcf>
#
# arguments:
#   ref Reference FASTA
#   in_vcf Input VCF (.vcf.gz)
#   out_vcf Output VCF (.vcf.gz, will be tabix-indexed)
#
# example:
#   bash scripts/fix_vcf.sh \
#   reference/hs37d5.fa \
#   archive/data_pipeline/19_all_hg002/synthetic_hg002.hg37.sv.vcf.gz \
#   archive/data_pipeline/19_all_hg002/synthetic_hg002.hg37.sv.fixed.vcf.gz
#
# dependencies: python3 (with pysam), bcftools, bgzip, tabix
set -euo pipefail

# logging

readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

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
in_vcf="$2"
out_vcf="$3"

[[ -f "$ref" ]]    || error "Reference not found: $ref"
[[ -f "$in_vcf" ]] || error "Input VCF not found: $in_vcf"

# dependency check

for cmd in python3 bcftools bgzip tabix; do
    command -v "$cmd" &>/dev/null || error "$cmd not found in PATH"
done

# main

info "Resolving symbolic alleles in $(basename "$in_vcf") ..."

zcat "$in_vcf" \
    | python3 "$SCRIPT_DIR/fix_vcf.py" "$ref" \
    | bcftools norm --check-ref s --fasta-ref "$ref" -N -m-any \
    | bgzip > "$out_vcf"

tabix "$out_vcf"

info "Done — $(basename "$out_vcf")"
