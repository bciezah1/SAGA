#!/bin/bash

if [ $# -lt 3 ]; then
  echo "Usage: $0 <input_path> <covariates> <qcovariates>"
  echo "Example: $0 ../../toy_data/ SEX,AGE,PC1,PC2,PC3 SEX"
  exit 1
fi

INPUT_PATH=$1
COVAR_LIST=$2
QCOVAR_LIST=$3

export PATH="$(pwd)/bin:$PATH"

# Call the analysis script
bash plink_pipeline.sh "$INPUT_PATH" "$COVAR_LIST" "$QCOVAR_LIST" 
