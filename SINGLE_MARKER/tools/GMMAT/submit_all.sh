#!/bin/bash

# Check arguments
if [ "$#" -ne 6 ]; then
  echo "Usage: $0 <working_dir> <KINSHIP_INPUT> <GENO_INPUT> <PHENO_INPUT> <MODEL> <TYPE>"
  exit 1
fi

WORKDIR="$1"
KINSHIP_INPUT="$2"
DOSAGE_INPUT="$3"
PHENO_INPUT="$4"
MODEL_FORMULA="$5"
TYPE=$6

mkdir -p "$WORKDIR/logs"

SECONDS=0

# Export variables so sub-scripts can use them
export WORKDIR KINSHIP_INPUT DOSAGE_INPUT PHENO_INPUT MODEL_FORMULA TYPE

# Start timer
START_TIME=$(date +%s)

echo "==> Step 1: Kinship + PCA"
bash step1_run_kinship_pca.sh

echo "==> Step 2: GWAS per chromosome"
bash step2_run_full_gwas_pipeline_noslurm.sh ${TYPE}

echo "==> Step 3: Postprocessing"
bash step3_postprocess_gwas_noslurm.sh

# Stop timer
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Print and log runtime
echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> runtime_log.txt


DURATION=$SECONDS
echo "? Step 2 completed in $(($DURATION / 60)) min $(($DURATION % 60)) sec"
