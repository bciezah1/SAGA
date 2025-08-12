#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <COVAR_LIST> <DIAGNOST_VAR> <TRAIT_TYPE>"
  exit 1
fi

GENO_INPUT=$1
PHENO_INPUT=$2
COVAR_LIST=$3
QCOVAR_LIST=$4
TRAIT_TYPE=$5

export PATH="$(pwd)/bin:$PATH"

# Start timer
START_TIME=$(date +%s)

# Run pipeline
bash plink_pipeline.sh "$GENO_INPUT" "$PHENO_INPUT" "$COVAR_LIST" "$QCOVAR_LIST" "$TRAIT_TYPE"

# Stop timer
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Print and log runtime
echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> runtime_log.txt

# organize
mkdir output
cd output
mkdir plots tables
cd ../
mv *jpg ./output/plots
mv manhattan_input.txt  pheno_with_pcs.txt  sum_stats.txt ./output/tables
rm mypc* asso*

