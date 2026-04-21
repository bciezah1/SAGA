#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
  cat <<EOF
Usage:
  ./run_pipeline_saige.sh <GENO_INPUT> <PHENO_FILE> <COVAR_LIST> <QCOVAR_LIST> <PHENO_COL> <TRAIT_TYPE> <WORKING_DIR>

Arguments:
  GENO_INPUT    PLINK prefix (without .bed/.bim/.fam)
  PHENO_FILE    Phenotype file with header, including FID and IID
  COVAR_LIST    Comma-separated categorical covariates for SAIGE covarColList
                Use NONE if not needed
  QCOVAR_LIST   Comma-separated quantitative covariates for SAIGE qCovarColList
                Use NONE if not needed
  PHENO_COL     Phenotype column to analyze
  TRAIT_TYPE    binary or quantitative
  WORKING_DIR   Output working directory

Example:
  ./run_pipeline_saige.sh /path/to/geno /path/to/pheno.txt SEX,COHORT AGE PHENO binary my_saige_run
EOF
}

if [ $# -ne 7 ]; then
  usage
  exit 1
fi

GENO_INPUT=$1
INPUT_PHENO=$2
COVAR_LIST=$3
QCOVAR_LIST=$4
PHENO_COL=$5
TRAIT_TYPE=$6
WORKING_DIR=$7

if [ ! -f "${GENO_INPUT}.bed" ] || [ ! -f "${GENO_INPUT}.bim" ] || [ ! -f "${GENO_INPUT}.fam" ]; then
  echo "ERROR: PLINK input files not found for prefix: ${GENO_INPUT}"
  exit 1
fi

if [ ! -f "${INPUT_PHENO}" ]; then
  echo "ERROR: Phenotype file not found: ${INPUT_PHENO}"
  exit 1
fi

mkdir -p "${WORKING_DIR}"

START_TIME=$(date +%s)

bash "${SCRIPT_DIR}/saige.sh" \
  "${GENO_INPUT}" \
  "${INPUT_PHENO}" \
  "${COVAR_LIST}" \
  "${QCOVAR_LIST}" \
  "${PHENO_COL}" \
  "${TRAIT_TYPE}" \
  "${WORKING_DIR}"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "Runtime: ${RUNTIME} seconds"
echo -e "DATE\tGENO_INPUT\tPHENO_FILE\tPHENO_COL\tTRAIT_TYPE\tRUNTIME_SECONDS" >> "${WORKING_DIR}/runtime_log.tsv"
echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t${GENO_INPUT}\t${INPUT_PHENO}\t${PHENO_COL}\t${TRAIT_TYPE}\t${RUNTIME}" >> "${WORKING_DIR}/runtime_log.tsv"