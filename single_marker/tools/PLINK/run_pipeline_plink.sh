#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 6 ]; then
  echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <COVAR_LIST> <DIAGNOST_VAR> <TRAIT_TYPE> <OUTPUT_DIR>"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GENO_INPUT="$1"
PHENO_INPUT="$2"
COVAR_LIST="$3"
DIAGNOST_VAR="$4"
TRAIT_TYPE="$5"
OUTPUT_DIR="$6"

export PLINK_BIN="${PLINK_BIN:-}"

START_TIME=$(date +%s)

bash "${SCRIPT_DIR}/plink_pipeline.sh" \
  "$GENO_INPUT" \
  "$PHENO_INPUT" \
  "$COVAR_LIST" \
  "$DIAGNOST_VAR" \
  "$TRAIT_TYPE" \
  "$OUTPUT_DIR"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> "${SCRIPT_DIR}/runtime_log.txt"