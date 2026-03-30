#!/usr/bin/env bash
# mendelian_concordance.sh — Mendelian concordance across preservation swaps
#
# merges single-sample DeepVariant VCFs into trios and runs bcftools
# mendelian2 to count Mendelian violations. Tests 5 trio combinations:
#
# 1. All Frozen (baseline)
# 2. All Ensilicated (full ensilicated control)
# 3. Swap Father (ensilicated father + frozen mother + frozen son)
# 4. Swap Mother (frozen father + ensilicated mother + frozen son)
# 5. Swap Son (frozen father + frozen mother + ensilicated son)
#
# usage:
#   bash scripts/mendelian_concordance.sh <vcf_dir> <out_dir>
#
# arguments:
#   vcf_dir Directory containing DeepVariant .vcf files
#   out_dir Output directory for merged VCFs and results
#
# example:
#   bash scripts/mendelian_concordance.sh \
#   archive/data_pipeline/05_variants \
#   results/mendelian
#
# dependencies: bcftools (with mendelian2 plugin), bgzip, tabix
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

vcf_dir="$1"
out_dir="$2"

[[ -d "$vcf_dir" ]] || error "VCF directory not found: $vcf_dir"

# dependency check

for cmd in bcftools bgzip tabix; do
    command -v "$cmd" &>/dev/null || error "$cmd not found in PATH"
done

# sample map (GIAB Ashkenazi trio)
#
# HG002 = Son (proband) Frozen: CD_3033 Ensilicated: CD_3032
# HG003 = Father Frozen: CD_3031 Ensilicated: CD_3030
# HG004 = Mother Frozen: CD_3029 Ensilicated: CD_3028

declare -A VCFS=(
    [HG002_frozen]="$vcf_dir/CD_3033_GIAB.nist.aligned.sorted.full.vcf"
    [HG002_ensilicated]="$vcf_dir/CD_3032_Cache.nist.aligned.sorted.full.vcf"
    [HG003_frozen]="$vcf_dir/CD_3031_GIAB.nist.aligned.sorted.full.vcf"
    [HG003_ensilicated]="$vcf_dir/CD_3030_Cache.nist.aligned.sorted.full.vcf"
    [HG004_frozen]="$vcf_dir/CD_3029_GIAB.nist.aligned.sorted.full.vcf"
    [HG004_ensilicated]="$vcf_dir/CD_3028_Cache.nist.aligned.sorted.full.vcf"
)

# prepare: rename samples, bgzip, index

mkdir -p "$out_dir/prep"

prepare_vcf() {
    local key="$1"
    local sample_name="$2"
    local src="${VCFS[$key]}"
    local dst="$out_dir/prep/${key}.vcf.gz"

    if [[ -f "$dst" && -f "${dst}.tbi" ]]; then
        info "  $key already prepared — skipping"
        return
    fi

    [[ -f "$src" ]] || error "VCF not found: $src"

    info "  Preparing $key (sample=$sample_name) ..."
    echo "$sample_name" > "$out_dir/prep/${key}.rename.txt"
    bcftools reheader -s "$out_dir/prep/${key}.rename.txt" "$src" \
        | bcftools view -Oz -o "$dst"
    tabix -p vcf "$dst"
}

info "Preparing individual VCFs ..."
for key in "${!VCFS[@]}"; do
    # extract sample name (HG002, HG003, HG004) from key
    sample_name="${key%%_*}"
    prepare_vcf "$key" "$sample_name"
done

# define trio combinations

declare -A TRIOS=(
    [all_frozen]="HG003_frozen,HG004_frozen,HG002_frozen"
    [all_ensilicated]="HG003_ensilicated,HG004_ensilicated,HG002_ensilicated"
    [swap_father]="HG003_ensilicated,HG004_frozen,HG002_frozen"
    [swap_mother]="HG003_frozen,HG004_ensilicated,HG002_frozen"
    [swap_son]="HG003_frozen,HG004_frozen,HG002_ensilicated"
)

# merge and run mendelian2 for each trio

mkdir -p "$out_dir/merged"

info ""
info "Running Mendelian concordance analysis ..."

results_file="$out_dir/mendelian_summary.tsv"
printf "trio\ttotal_sites\tconsistent\tinconsistent\tmissing\tconcordance_rate\n" > "$results_file"

for trio_name in all_frozen all_ensilicated swap_father swap_mother swap_son; do
    IFS=',' read -r father mother son <<< "${TRIOS[$trio_name]}"

    father_vcf="$out_dir/prep/${father}.vcf.gz"
    mother_vcf="$out_dir/prep/${mother}.vcf.gz"
    son_vcf="$out_dir/prep/${son}.vcf.gz"

    merged_vcf="$out_dir/merged/${trio_name}.vcf.gz"
    mendel_out="$out_dir/${trio_name}.mendelian.txt"

    info ""
    info "--- $trio_name ---"
    info "  Father: $father"
    info "  Mother: $mother"
    info "  Son:    $son"

    # merge three single-sample VCFs into one multi-sample VCF
    if [[ ! -f "$merged_vcf" ]]; then
        info "  Merging ..."
        bcftools merge --force-samples \
            "$father_vcf" "$mother_vcf" "$son_vcf" \
            -Oz -o "$merged_vcf"
        tabix -p vcf "$merged_vcf"
    fi

    # run mendelian2 (HG002 is male proband, hence 1X:)
    info "  Running mendelian2 ..."
    bcftools +mendelian2 "$merged_vcf" \
        -p "1X:HG002,HG003,HG004" \
        -m c \
        > "$mendel_out" 2>&1

    info "  Results:"
    cat "$mendel_out" | while IFS= read -r line; do info "    $line"; done

    # parse counts from mendelian2 tab-delimited output
    consistent=$(grep -P '^ngood\t' "$mendel_out" | cut -f2 || echo "0")
    inconsistent=$(grep -P '^nmerr\t' "$mendel_out" | cut -f2 || echo "0")
    missing=$(grep -P '^nmissing\t' "$mendel_out" | cut -f2 || echo "0")
    consistent=${consistent:-0}
    inconsistent=${inconsistent:-0}
    missing=${missing:-0}
    total=$((consistent + inconsistent))

    if [[ "$total" -gt 0 ]]; then
        rate=$(awk "BEGIN {printf \"%.6f\", $consistent / $total}")
    else
        rate="NA"
    fi

    printf "%s\t%d\t%d\t%d\t%d\t%s\n" \
        "$trio_name" "$total" "$consistent" "$inconsistent" "$missing" "$rate" \
        >> "$results_file"
done

info ""

info "Summary:"
info ""
column -t -s$'\t' "$results_file"
info ""
info "Full results: $results_file"
info "Done."
