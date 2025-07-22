#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Usage: $0 <input_path> <covariates> <qcovariates> <model_tag> <chr>"
  echo "Example: $0 ../../toy_data/ SEX,AGE,PC1,PC2,PC3 SEX model1 1"
  exit 1
fi

INPUT_PATH=$1
COVAR_LIST=$2
QCOVAR_LIST=$3
MODEL_TAG=$4
CHR=$5  # Single chromosome to process (e.g. 1 to 22)

export PATH="$(pwd)/bin:$PATH"

# Call the analysis script
bash saige.sh "$INPUT_PATH" "$COVAR_LIST" "$QCOVAR_LIST" "$MODEL_TAG" "$CHR"

# how to run it:
#  ./run_pipeline_saige.sh /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/PLINK_FILES/PRADI/ AGE,SEX SEX model1 1
