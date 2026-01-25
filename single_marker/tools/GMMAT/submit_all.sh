#!/bin/bash
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <working_dir> <GENO_INPUT> <PHENO_INPUT> <MODEL> <TYPE> <OUTPUT_DIR>"
    exit 1
fi

#WORKDIR="$1"
GENO_INPUT="$1"
PHENO_INPUT="$2"
MODEL_FORMULA="$3"
TYPE="$4"
OUTPUT_DIR="$5"

mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/tmp"

export GENO_INPUT PHENO_INPUT MODEL_FORMULA TYPE OUTPUT_DIR

SECONDS=0
START_TIME=$(date +%s)

echo "==> Step 1: Kinship + PCA"
bash step1_run_kinship_pca.sh

echo "==> Step 2: GWAS per chromosome"
bash step2_run_full_gwas_pipeline_noslurm.sh "$TYPE"

echo "==> Step 3: Postprocessing"
bash step3_postprocess_gwas_noslurm.sh

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> "$OUTPUT_DIR/runtime_log.txt"
cd "$OUTPUT_DIR"
rm -r assoc_step* mypc.* QCed.* pruned_snps.* kinship_step1.* mykinship.* *.tmp mega_scores_chr1_kinship_09_23_2024_5PCs_sorted_model1.txt pheno_f* tmp logs 
mkdir output
cd output
mkdir plots tables
cd ../
mv manhattan_input.txt pheno_with_pcs.txt sum_stats.txt ./output/tables
cd ../
mv *jpg "$OUTPUT_DIR/output/plots"

rm -r output
