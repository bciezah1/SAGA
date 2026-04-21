#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PLINK_BIN="${SCRIPT_DIR}/../bin/plink"
DOWNSAMPLE_R="${SCRIPT_DIR}/step4_downsampling.R"

ASSOCIATION_MAF_THRESHOLD=${ASSOCIATION_MAF_THRESHOLD:-0.001}
CLEAN_TEMP_QC=${CLEAN_TEMP_QC:-false}

usage() {
  cat <<EOF
Usage:
  bash qc.sh <input_gwas_prefix> <input_kinship_prefix> [assoc_output_prefix] [kinship_output_prefix]

Example:
  bash qc.sh /path/to/geno /path/to/geno QCed.assoc QCed.kinship

Environment variables:
  ASSOCIATION_MAF_THRESHOLD   MAF threshold for association data (default: 0.001)
  CLEAN_TEMP_QC               true or false; remove temporary QC files at end (default: false)
EOF
}

if [ $# -lt 2 ] || [ $# -gt 4 ]; then
  usage
  exit 1
fi

INPUT_GWAS=$1
INPUT_KINSHIP=$2
ASSOC_OUT=${3:-QCed.assoc}
KINSHIP_OUT=${4:-QCed.kinship}

for prefix in "${INPUT_GWAS}" "${INPUT_KINSHIP}"; do
  if [ ! -f "${prefix}.bed" ] || [ ! -f "${prefix}.bim" ] || [ ! -f "${prefix}.fam" ]; then
    echo "ERROR: Missing PLINK files for prefix: ${prefix}"
    exit 1
  fi
done

if [ ! -x "${PLINK_BIN}" ]; then
  echo "ERROR: PLINK binary not found or not executable: ${PLINK_BIN}"
  exit 1
fi

if [ ! -f "${DOWNSAMPLE_R}" ]; then
  echo "ERROR: Downsampling R script not found: ${DOWNSAMPLE_R}"
  exit 1
fi

have_plink_prefix() {
  local prefix=$1
  [ -s "${prefix}.bed" ] && [ -s "${prefix}.bim" ] && [ -s "${prefix}.fam" ]
}

count_samples() {
  wc -l < "${1}.fam"
}

count_variants() {
  wc -l < "${1}.bim"
}

append_summary() {
  local step=$1
  local dataset=$2
  local prefix=$3
  echo -e "${step}\t${dataset}\t$(count_samples "${prefix}")\t$(count_variants "${prefix}")" >> "${QC_SUMMARY}"
}

QC_SUMMARY="qc_summary.tsv"
echo -e "STEP\tDATASET\tN_SAMPLES\tN_VARIANTS" > "${QC_SUMMARY}"

echo "Starting QC pipeline..."
echo "GWAS input prefix           : ${INPUT_GWAS}"
echo "Kinship input prefix        : ${INPUT_KINSHIP}"
echo "Association MAF threshold   : ${ASSOCIATION_MAF_THRESHOLD}"
echo "Clean temporary QC files    : ${CLEAN_TEMP_QC}"

append_summary "INPUT" "GWAS" "${INPUT_GWAS}"
append_summary "INPUT" "KINSHIP" "${INPUT_KINSHIP}"

echo "Step 1. SNP missingness filter for association data..."
if have_plink_prefix "assoc_step1"; then
  echo "  Found existing assoc_step1.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile "${INPUT_GWAS}" \
    --geno 0.02 \
    --make-bed \
    --out assoc_step1
fi
append_summary "STEP1_GENO" "assoc_step1" "assoc_step1"

echo "Step 2. Individual missingness filter for association data..."
if have_plink_prefix "assoc_step2"; then
  echo "  Found existing assoc_step2.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile assoc_step1 \
    --mind 0.02 \
    --make-bed \
    --out assoc_step2
fi
append_summary "STEP2_MIND" "assoc_step2" "assoc_step2"

echo "Step 3. Minor allele frequency filter for association data..."
if have_plink_prefix "assoc_step3"; then
  echo "  Found existing assoc_step3.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile assoc_step2 \
    --maf "${ASSOCIATION_MAF_THRESHOLD}" \
    --make-bed \
    --out assoc_step3
fi
append_summary "STEP3_MAF" "assoc_step3" "assoc_step3"

echo "Step 4. Hardy-Weinberg equilibrium filter for association data..."
if have_plink_prefix "${ASSOC_OUT}"; then
  echo "  Found existing ${ASSOC_OUT}.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile assoc_step3 \
    --hwe 1e-6 \
    --make-bed \
    --out "${ASSOC_OUT}"
fi
append_summary "STEP4_HWE" "${ASSOC_OUT}" "${ASSOC_OUT}"

echo "Association QC complete. Output: ${ASSOC_OUT}.*"

echo "Step 5. LD pruning for kinship data..."
if [ -s "temp_snps.prune.in" ]; then
  echo "  Found existing temp_snps.prune.in ; skipping."
else
  "${PLINK_BIN}" \
    --bfile "${INPUT_KINSHIP}" \
    --indep-pairwise 50 5 0.2 \
    --out temp_snps
fi

echo "Step 6. Extract pruned SNPs for kinship data..."
if have_plink_prefix "temp_snps_dataset"; then
  echo "  Found existing temp_snps_dataset.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile "${INPUT_KINSHIP}" \
    --extract temp_snps.prune.in \
    --make-bed \
    --out temp_snps_dataset
fi
append_summary "STEP6_LD_PRUNED" "temp_snps_dataset" "temp_snps_dataset"

echo "Step 7. Compute allele frequencies for downsampling..."
if [ -s "temp_snps_dataset.freq.frq" ]; then
  echo "  Found existing temp_snps_dataset.freq.frq ; skipping."
else
  "${PLINK_BIN}" \
    --bfile temp_snps_dataset \
    --freq \
    --out temp_snps_dataset.freq
fi

echo "Step 8. Downsample SNPs by MAF bin for sparse GRM..."
if [ -s "list_snps_for_grm.txt" ]; then
  echo "  Found existing list_snps_for_grm.txt ; skipping."
else
  Rscript "${DOWNSAMPLE_R}" \
    "temp_snps_dataset.freq.frq" \
    "list_snps_for_grm.txt"
fi

if [ ! -f "list_snps_for_grm.txt" ]; then
  echo "ERROR: list_snps_for_grm.txt was not generated."
  exit 1
fi

if [ ! -s "list_snps_for_grm.txt" ]; then
  echo "ERROR: list_snps_for_grm.txt is empty."
  exit 1
fi

echo "Step 9. Extract downsampled SNPs for kinship data..."
if have_plink_prefix "${KINSHIP_OUT}"; then
  echo "  Found existing ${KINSHIP_OUT}.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile temp_snps_dataset \
    --extract list_snps_for_grm.txt \
    --make-bed \
    --out "${KINSHIP_OUT}"
fi
append_summary "STEP9_DOWNSAMPLED" "${KINSHIP_OUT}" "${KINSHIP_OUT}"

echo "Kinship QC complete. Output: ${KINSHIP_OUT}.*"
echo "QC summary written to: ${QC_SUMMARY}"

if [ "${CLEAN_TEMP_QC}" = "true" ]; then
  echo "Cleaning temporary QC files..."
  rm -f assoc_step1.bed assoc_step1.bim assoc_step1.fam assoc_step1.log assoc_step1.nosex
  rm -f assoc_step2.bed assoc_step2.bim assoc_step2.fam assoc_step2.log assoc_step2.nosex
  rm -f assoc_step3.bed assoc_step3.bim assoc_step3.fam assoc_step3.log assoc_step3.nosex
  rm -f temp_snps.prune.in temp_snps.prune.out temp_snps.log temp_snps.nosex
  rm -f temp_snps_dataset.bed temp_snps_dataset.bim temp_snps_dataset.fam temp_snps_dataset.log temp_snps_dataset.nosex
  rm -f temp_snps_dataset.freq.frq temp_snps_dataset.freq.log temp_snps_dataset.freq.nosex
else
  echo "Keeping temporary QC files to support resume behavior."
fi

echo "QC pipeline completed successfully."
echo "Final outputs:"
echo "  ${ASSOC_OUT}.*"
echo "  ${KINSHIP_OUT}.*"
echo "  list_snps_for_grm.txt"
echo "  ${QC_SUMMARY}"