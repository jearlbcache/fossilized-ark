#!/usr/bin/env bash
# per_read_identity.sh — compute per-read identity histograms from BAMs
#
# uses samtools view + awk to extract NM tag and aligned length from CIGAR,
# bins identity into 0.01% bins on-the-fly. outputs a histogram TSV per BAM.
#
# usage:
#   bash scripts/per_read_identity.sh
#
# dependencies: samtools
set -euo pipefail

export PATH="/home/anaconda3/bin:$PATH"

PROJECT="/data/jb/project/giab_stanford"
BAM_DIR="$PROJECT/archive/data_pipeline/12_giab_hg38_aligned"
OUT_DIR="$PROJECT/results/error_profiles"
mkdir -p "$OUT_DIR"

declare -A BAMS=(
    [HG002_frozen]="$BAM_DIR/CD_3033_GIAB.nist.aligned.sorted.bam"
    [HG002_ensilicated]="$BAM_DIR/CD_3032_Cache.nist.aligned.sorted.bam"
    [HG003_frozen]="$BAM_DIR/CD_3031_GIAB.nist.aligned.sorted.bam"
    [HG003_ensilicated]="$BAM_DIR/CD_3030_Cache.nist.aligned.sorted.bam"
    [HG004_frozen]="$BAM_DIR/CD_3029_GIAB.nist.aligned.sorted.bam"
    [HG004_ensilicated]="$BAM_DIR/CD_3028_Cache.nist.aligned.sorted.bam"
)

run_one() {
    local label="$1"
    local bam="${BAMS[$label]}"
    local out="$OUT_DIR/${label}_identity_hist.tsv"
    echo "$(date +%H:%M:%S) START $label"

    samtools view -@ 14 -F 0x904 "$bam" | awk -F'\t' '
    {
        nm=-1
        for(i=12;i<=NF;i++) {
            if(substr($i,1,5)=="NM:i:") {
                nm=substr($i,6)+0
                break
            }
        }
        if(nm<0) next

        cigar=$6; alen=0
        while(match(cigar,/[0-9]+/,a)) {
            n=a[0]+0
            op=substr(cigar,RSTART+RLENGTH,1)
            if(op=="M"||op=="="||op=="X"||op=="I"||op=="D") alen+=n
            cigar=substr(cigar,RSTART+RLENGTH+1)
        }
        if(alen>0) {
            bin=int((1-nm/alen)*10000)
            if(bin<0) bin=0
            hist[bin]++
        }
    }
    END {
        for(bin in hist) printf "%d\t%d\n", bin, hist[bin]
    }' | sort -n > "$out"

    echo "$(date +%H:%M:%S) DONE  $label"
}

# run all 6 in parallel (14 threads each = 84 of 96 threads)
for label in HG002_frozen HG002_ensilicated HG003_frozen HG003_ensilicated HG004_frozen HG004_ensilicated; do
    run_one "$label" &
done
wait

echo "all done."
