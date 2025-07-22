#!/bin/bash

# Check arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <working_dir> <BFILE_PATH_MAF0_05_R2_099> <DOSAGE_DIR_MAF0_01_R2_080> <MODEL_TYPE>"
  exit 1
fi

WORKDIR="$1"
BFILE_PATH_MAF0_05_R2_099="$2"
DOSAGE_DIR_MAF0_01_R2_080="$3"
MODEL_FORMULA="$4"

mkdir -p "$WORKDIR/logs"

SECONDS=0

# Export variables so sub-scripts can use them
export WORKDIR BFILE_PATH_MAF0_05_R2_099 DOSAGE_DIR_MAF0_01_R2_080 MODEL_FORMULA

echo "==> Step 1: Kinship + PCA"
bash step1_run_kinship_pca.sh

echo "==> Step 2: GWAS per chromosome"
bash step2_run_full_gwas_pipeline_noslurm.sh

echo "==> Step 3: Postprocessing"
bash step3_postprocess_gwas_noslurm.sh

DURATION=$SECONDS
echo "? Step 2 completed in $(($DURATION / 60)) min $(($DURATION % 60)) sec"
