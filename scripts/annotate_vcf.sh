#!/usr/bin/env bash
# annotate_vcf.sh — Resolve symbolic SV alleles in a VCF
#
# wrapper around fix_vcf.sh — resolves <DEL> alleles with the reference
# sequence, normalizes, compresses, and indexes the output.
#
# usage:
#   bash scripts/annotate_vcf.sh <ref> <input_vcf> <output_vcf>
#
# arguments:
#   ref Reference FASTA
#   input_vcf Input SV VCF (.vcf.gz)
#   output_vcf Output annotated VCF (.vcf.gz)
#
# example:
#   bash scripts/annotate_vcf.sh \
#   reference/hs37d5.fa \
#   archive/data_pipeline/19_all_hg002/synthetic_hg002.hg37.sv.vcf.gz \
#   archive/data_pipeline/19_all_hg002/synthetic_hg002.hg37.sv.annotated.vcf.gz
#
# dependencies: fix_vcf.sh (same directory)
set -euo pipefail

# logging

readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

info()  { printf '[%s] %s\n' "$SCRIPT_NAME" "$*"; }
error() { printf '[%s] ERROR: %s\n' "$SCRIPT_NAME" "$*" >&2; exit 1; }

# arguments

usage() {
    sed -n '/^# Usage:/,/^# Dependencies:/{ /^# Dependencies:/d; s/^# //; s/^#//; p }' "$0"
    exit 1
}

[[ $# -ge 3 ]] || usage

ref="$1"
input_vcf="$2"
output_vcf="$3"

[[ -f "$ref" ]]       || error "Reference not found: $ref"
[[ -f "$input_vcf" ]] || error "Input VCF not found: $input_vcf"

# main

bash "$SCRIPT_DIR/fix_vcf.sh" "$ref" "$input_vcf" "$output_vcf"
