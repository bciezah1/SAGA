#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SAGA_GMMAT_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

model="model1"

RAW_GWAS_FILE="$OUTPUT_DIR/mega_scores_chr1_${model}_raw.txt"
MANHATTAN_INPUT="$OUTPUT_DIR/manhattan_input.txt"
SUM_STATS="$OUTPUT_DIR/sum_stats.txt"
PLOT_DIR="$OUTPUT_DIR/output/plots"

TEMP1="$OUTPUT_DIR/temp1.tmp"
TEMP_FULL="$OUTPUT_DIR/temp_full.tmp"
TEMP_FILTERED="$OUTPUT_DIR/temp_filtered.tmp"

mkdir -p "$PLOT_DIR"

if [ ! -f "$RAW_GWAS_FILE" ]; then
    echo "ERROR: expected GWAS result file not found: $RAW_GWAS_FILE"
    exit 1
fi

awk 'NR > 1 {print $2 "\t" $1 "\t" $4 "\t" $11}' "$RAW_GWAS_FILE" > "$TEMP1"

echo -e "SNP\tCHR\tBP\tPVAL" > "$MANHATTAN_INPUT"
cat "$TEMP1" >> "$MANHATTAN_INPUT"

awk -v OFS="\t" '{
    beta = ($10 + 0 != 0) ? $9 / $10 : "NA"
    se   = ($10 + 0 != 0) ? 1 / sqrt($10) : "NA"
    split($2, loc, ":")
    print loc[1], loc[2], $0, beta, se
}' "$RAW_GWAS_FILE" > "$TEMP_FULL"

awk -F'\t' '$10 >= 0.01 && $10 < 0.99' "$TEMP_FULL" > "$TEMP_FILTERED"

echo -e "CHR\tBP\tCHR_ORIG\tSNP\tCM\tBP_ORIG\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\tBETA\tSE" > "$SUM_STATS"
cat "$TEMP_FILTERED" >> "$SUM_STATS"

Rscript "$SCRIPT_DIR/create_manhattan.R" "$MANHATTAN_INPUT" "$PLOT_DIR"
Rscript "$SCRIPT_DIR/create_qq.plot.R" "$MANHATTAN_INPUT" "$PLOT_DIR"
Rscript "$SCRIPT_DIR/create_circular_manhattan.R" "$MANHATTAN_INPUT" "$PLOT_DIR"
Rscript "$SCRIPT_DIR/create_density_plot.R" "$MANHATTAN_INPUT" "$PLOT_DIR"

rm -f "$TEMP1" "$TEMP_FULL" "$TEMP_FILTERED"

echo "Post-processing completed in $OUTPUT_DIR"