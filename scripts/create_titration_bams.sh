#!/usr/bin/env bash
# create_titration_bams.sh — Create ensilicated/frozen titration BAMs for HG002
#
# subsamples and merges ensilicated (CD_3032) and frozen (CD_3033) BAMs at 7
# mixture ratios: 0%, 10%, 25%, 50%, 75%, 90%, 100% ensilicated.
#
# usage:
#   bash scripts/create_titration_bams.sh [--only PCT] [--dry-run]
#
# options:
#   --only PCT Create only the specified percentage (e.g., --only 10)
#   --dry-run Print commands without executing
#
# dependencies: samtools, python3
set -euo pipefail

export PATH="/home/anaconda3/bin:$PATH"

readonly SCRIPT_NAME="$(basename "$0")"
info()  { printf '[%s] %s\n' "$SCRIPT_NAME" "$*"; }
error() { printf '[%s] ERROR: %s\n' "$SCRIPT_NAME" "$*" >&2; exit 1; }

# configuration

PROJECT_ROOT="/data/jb/project/giab_stanford"
BAM_DIR="$PROJECT_ROOT/archive/data_pipeline/12_giab_hg38_aligned"
OUT="$PROJECT_ROOT/results/titration"

FROZEN_BAM="$BAM_DIR/CD_3033_GIAB.nist.aligned.sorted.bam"
ENSILICATED_BAM="$BAM_DIR/CD_3032_Cache.nist.aligned.sorted.bam"

FROZEN_TOTAL=33346590
ENSILICATED_TOTAL=24098474
TARGET=12000000

THREADS=16

# parse arguments

ONLY=""
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only)   ONLY="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        *)        error "Unknown argument: $1" ;;
    esac
done

# helpers

# compute samtools -s SEED.FRAC string.
# samtools interprets -s 42.1234 as seed=42, keep fraction=0.1234
# critical: use [42 + fraction] so the decimal digits ARE the fraction.
# the old code used f'42.{f*10000:.0f}' which dropped leading zeros
# (e.g., 0.0694 -> 42.694 = 69.4% instead of 42.0694 = 6.94%).
compute_subsample_arg() {
    local n_reads="$1"
    local total="$2"
    python3 -c "print(f'{42 + $n_reads/$total:.4f}')"
}

run_cmd() {
    if $DRY_RUN; then
        info "[dry-run] $*"
    else
        "$@"
    fi
}

# validation

[[ -f "$FROZEN_BAM" ]] || error "Frozen BAM not found: $FROZEN_BAM"
[[ -f "$ENSILICATED_BAM" ]] || error "Ensilicated BAM not found: $ENSILICATED_BAM"
mkdir -p "$OUT"

# create titration BAMs

ALL_PCTS=(0 10 25 50 75 90 100)

if [[ -n "$ONLY" ]]; then
    PCTS=("$ONLY")
else
    PCTS=("${ALL_PCTS[@]}")
fi

for pct_ensilicated in "${PCTS[@]}"; do
    pct_frozen=$((100 - pct_ensilicated))
    n_ensilicated=$((TARGET * pct_ensilicated / 100))
    n_frozen=$((TARGET * pct_frozen / 100))

    out_bam="$OUT/hg002_mix_${pct_ensilicated}pct_ensilicated.bam"

    sfrozen=$(compute_subsample_arg "$n_frozen" "$FROZEN_TOTAL")
    sensilicated=$(compute_subsample_arg "$n_ensilicated" "$ENSILICATED_TOTAL")

    info "=== ${pct_ensilicated}% ensilicated (${pct_frozen}% frozen): sfrozen=${sfrozen} sensilicated=${sensilicated} ==="

    if [[ "$pct_ensilicated" -eq 0 ]]; then
        run_cmd samtools view -@ "$THREADS" -bs "$sfrozen" "$FROZEN_BAM" \
            | run_cmd samtools sort -@ "$THREADS" -o "$out_bam"
    elif [[ "$pct_ensilicated" -eq 100 ]]; then
        run_cmd samtools view -@ "$THREADS" -bs "$sensilicated" "$ENSILICATED_BAM" \
            | run_cmd samtools sort -@ "$THREADS" -o "$out_bam"
    else
        run_cmd samtools view -@ 8 -bs "$sfrozen" "$FROZEN_BAM" -o "$OUT/tmp_frozen.bam" &
        run_cmd samtools view -@ 8 -bs "$sensilicated" "$ENSILICATED_BAM" -o "$OUT/tmp_ensilicated.bam" &
        wait
        run_cmd samtools merge -@ "$THREADS" -f "$OUT/tmp_merged.bam" "$OUT/tmp_frozen.bam" "$OUT/tmp_ensilicated.bam"
        run_cmd samtools sort -@ "$THREADS" -o "$out_bam" "$OUT/tmp_merged.bam"
        rm -f "$OUT/tmp_frozen.bam" "$OUT/tmp_ensilicated.bam" "$OUT/tmp_merged.bam"
    fi

    run_cmd samtools index -@ "$THREADS" "$out_bam"

    if ! $DRY_RUN; then
        reads=$(samtools view -c "$out_bam")
        info "  Done: ${reads} reads (target: ~$((n_frozen + n_ensilicated)))"

        # sanity check: reads should be within 15% of target
        expected=$((n_frozen + n_ensilicated))
        lower=$((expected * 85 / 100))
        upper=$((expected * 115 / 100))
        if [[ "$reads" -lt "$lower" || "$reads" -gt "$upper" ]]; then
            error "Read count $reads is outside ±15% of expected $expected for ${pct_ensilicated}% ensilicated!"
        fi
    fi
    info ""
done

info "All requested titration BAMs created."
