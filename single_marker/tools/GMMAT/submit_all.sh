#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <MODEL> <TYPE> <OUTPUT_DIR>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GENO_INPUT="$1"
PHENO_INPUT="$2"
MODEL_FORMULA="$3"
TYPE="$4"
OUTPUT_DIR="$5"

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/tmp"

export GENO_INPUT
export PHENO_INPUT
export MODEL_FORMULA
export TYPE
export OUTPUT_DIR
export SAGA_GMMAT_SCRIPT_DIR="$SCRIPT_DIR"

START_TIME=$(date +%s)

echo "==> Step 1: Kinship + PCA"
bash "$SCRIPT_DIR/step1_run_kinship_pca.sh"

echo "==> Step 2: GWAS per chromosome"
bash "$SCRIPT_DIR/step2_run_full_gwas_pipeline_noslurm.sh" "$TYPE"

echo "==> Step 3: Postprocessing"
bash "$SCRIPT_DIR/step3_postprocess_gwas_noslurm.sh"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> "$OUTPUT_DIR/runtime_log.txt"

mkdir -p "$OUTPUT_DIR/output/plots"
mkdir -p "$OUTPUT_DIR/output/tables"

for f in manhattan_input.txt pheno_with_pcs.txt sum_stats.txt; do
    if [ -f "$OUTPUT_DIR/$f" ]; then
        mv "$OUTPUT_DIR/$f" "$OUTPUT_DIR/output/tables/"
    fi
done

shopt -s nullglob
jpg_files=("$OUTPUT_DIR"/*.jpg)
if [ "${#jpg_files[@]}" -gt 0 ]; then
    mv "${jpg_files[@]}" "$OUTPUT_DIR/output/plots/"
fi
shopt -u nullglob

rm -f "$OUTPUT_DIR"/assoc_step1.bed "$OUTPUT_DIR"/assoc_step1.bim "$OUTPUT_DIR"/assoc_step1.fam
rm -f "$OUTPUT_DIR"/assoc_step1.log "$OUTPUT_DIR"/assoc_step1.nosex
rm -f "$OUTPUT_DIR"/assoc_step2.bed "$OUTPUT_DIR"/assoc_step2.bim "$OUTPUT_DIR"/assoc_step2.fam
rm -f "$OUTPUT_DIR"/assoc_step2.log "$OUTPUT_DIR"/assoc_step2.nosex
rm -f "$OUTPUT_DIR"/assoc_step3.bed "$OUTPUT_DIR"/assoc_step3.bim "$OUTPUT_DIR"/assoc_step3.fam
rm -f "$OUTPUT_DIR"/assoc_step3.log "$OUTPUT_DIR"/assoc_step3.nosex

rm -f "$OUTPUT_DIR"/kinship_step1.bed "$OUTPUT_DIR"/kinship_step1.bim "$OUTPUT_DIR"/kinship_step1.fam
rm -f "$OUTPUT_DIR"/kinship_step1.log "$OUTPUT_DIR"/kinship_step1.nosex

rm -f "$OUTPUT_DIR"/pruned_snps.prune.in "$OUTPUT_DIR"/pruned_snps.prune.out "$OUTPUT_DIR"/pruned_snps.log

rm -f "$OUTPUT_DIR"/*.tmp
rm -f "$OUTPUT_DIR"/pheno_for_kinship.txt
rm -f "$OUTPUT_DIR"/pheno_fid_kinship.txt

rm -rf "$OUTPUT_DIR/tmp"
rm -rf "$OUTPUT_DIR/logs"