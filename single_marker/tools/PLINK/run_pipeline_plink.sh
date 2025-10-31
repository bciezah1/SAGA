#!/bin/bash

if [ $# -lt 6 ]; then
  echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <COVAR_LIST> <DIAGNOST_VAR> <TRAIT_TYPE> <OUTPUT_DIR>"
  exit 1
fi

GENO_INPUT=$1
PHENO_INPUT=$2
COVAR_LIST=$3
QCOVAR_LIST=$4
TRAIT_TYPE=$5
OUTPUT_DIR=$6       

export PATH="$(pwd)/bin:$PATH"

# Start timer
START_TIME=$(date +%s)

# Run pipeline
bash plink_pipeline.sh "$GENO_INPUT" "$PHENO_INPUT" "$COVAR_LIST" "$QCOVAR_LIST" "$TRAIT_TYPE" "$OUTPUT_DIR"

# Stop timer
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Print and log runtime
echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> runtime_log.txt


