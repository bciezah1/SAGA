#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 2 ]; then
  echo "Usage: $0 <GENO_INPUT_FOR_ASSOC> <GENO_INPUT_FOR_KINSHIP>"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_GWAS="$1"
INPUT_KINSHIP="$2"

ASSOCIATION_MAF_THRESHOLD="${ASSOCIATION_MAF_THRESHOLD:-0.01}"

resolve_plink() {
  if [ -n "${PLINK_BIN:-}" ] && [ -x "${PLINK_BIN}" ]; then
    echo "${PLINK_BIN}"
    return 0
  fi

  if command -v plink >/dev/null 2>&1; then
    command -v plink
    return 0
  fi

  if [ -x "${SCRIPT_DIR}/../bin/plink" ]; then
    echo "${SCRIPT_DIR}/../bin/plink"
    return 0
  fi

  echo "ERROR: Could not find PLINK. Set PLINK_BIN or add plink to PATH." >&2
  exit 1
}

PLINK="$(resolve_plink)"

plink_prefix_exists() {
  local prefix="$1"
  [ -f "${prefix}.bed" ] && [ -f "${prefix}.bim" ] && [ -f "${prefix}.fam" ]
}

require_input_prefix() {
  local prefix="$1"
  if [ ! -f "${prefix}.bed" ] || [ ! -f "${prefix}.bim" ] || [ ! -f "${prefix}.fam" ]; then
    echo "ERROR: Missing PLINK files for prefix: ${prefix}" >&2
    exit 1
  fi
}

count_samples() {
  wc -l < "${1}.fam"
}

count_variants() {
  wc -l < "${1}.bim"
}

append_summary() {
  local step="$1"
  local dataset="$2"
  local prefix="$3"
  if plink_prefix_exists "${prefix}"; then
    echo -e "${step}\t${dataset}\t$(count_samples "${prefix}")\t$(count_variants "${prefix}")" >> "${QC_SUMMARY}"
  fi
}

require_input_prefix "${INPUT_GWAS}"
require_input_prefix "${INPUT_KINSHIP}"

QC_SUMMARY="qc_summary.tsv"
echo -e "STEP\tDATASET\tN_SAMPLES\tN_VARIANTS" > "${QC_SUMMARY}"

echo "Using PLINK: ${PLINK}"
echo "GWAS input prefix         : ${INPUT_GWAS}"
echo "Kinship input prefix      : ${INPUT_KINSHIP}"
echo "Association MAF threshold : ${ASSOCIATION_MAF_THRESHOLD}"

append_summary "INPUT" "GWAS" "${INPUT_GWAS}"
append_summary "INPUT" "KINSHIP" "${INPUT_KINSHIP}"

# ==========================================
# Association QC
# ==========================================

echo "Step 1. SNP missingness filter for association data..."
if plink_prefix_exists "assoc_step1"; then
  echo "  Found existing assoc_step1.* ; skipping."
else
  "${PLINK}" \
    --bfile "${INPUT_GWAS}" \
    --geno 0.02 \
    --make-bed \
    --out assoc_step1
fi
append_summary "STEP1_GENO" "assoc_step1" "assoc_step1"

echo "Step 2. Individual missingness filter for association data..."
if plink_prefix_exists "assoc_step2"; then
  echo "  Found existing assoc_step2.* ; skipping."
else
  "${PLINK}" \
    --bfile assoc_step1 \
    --mind 0.02 \
    --make-bed \
    --out assoc_step2
fi
append_summary "STEP2_MIND" "assoc_step2" "assoc_step2"

echo "Step 3. Minor allele frequency filter for association data..."
if plink_prefix_exists "assoc_step3"; then
  echo "  Found existing assoc_step3.* ; skipping."
else
  "${PLINK}" \
    --bfile assoc_step2 \
    --maf "${ASSOCIATION_MAF_THRESHOLD}" \
    --make-bed \
    --out assoc_step3
fi
append_summary "STEP3_MAF" "assoc_step3" "assoc_step3"

echo "Step 4. Hardy-Weinberg equilibrium filter for association data..."
if plink_prefix_exists "QCed.assoc"; then
  echo "  Found existing QCed.assoc.* ; skipping."
else
  "${PLINK}" \
    --bfile assoc_step3 \
    --hwe 1e-6 \
    --make-bed \
    --out QCed.assoc
fi
append_summary "STEP4_HWE" "QCed.assoc" "QCed.assoc"

echo "Association QC complete. Output: QCed.assoc.*"

# ==========================================
# PCA input QC (QCed.kinship)
# Use GMMAT-style parameters, but do not compute kinship matrix
# ==========================================

echo "Step 5. SNP/sample QC for PCA input (QCed.kinship)..."
if plink_prefix_exists "kinship_step1"; then
  echo "  Found existing kinship_step1.* ; skipping."
else
  "${PLINK}" \
    --bfile "${INPUT_KINSHIP}" \
    --geno 0.05 \
    --maf 0.05 \
    --mind 0.15 \
    --hwe 1e-6 \
    --make-bed \
    --out kinship_step1
fi
append_summary "STEP5_PCA_QC" "kinship_step1" "kinship_step1"

echo "Step 6. LD pruning for PCA input..."
if [ -f "pruned_snps.prune.in" ] && [ -f "pruned_snps.prune.out" ]; then
  echo "  Found existing pruned_snps.prune.in / pruned_snps.prune.out ; skipping."
else
  "${PLINK}" \
    --bfile kinship_step1 \
    --indep-pairwise 50 5 0.2 \
    --out pruned_snps
fi

if [ ! -f "pruned_snps.prune.in" ]; then
  echo "ERROR: pruned_snps.prune.in was not generated." >&2
  exit 1
fi

if [ ! -s "pruned_snps.prune.in" ]; then
  echo "ERROR: pruned_snps.prune.in is empty." >&2
  exit 1
fi

echo "Step 7. Extract LD-pruned SNPs for PCA input..."
if plink_prefix_exists "QCed.kinship"; then
  echo "  Found existing QCed.kinship.* ; skipping."
else
  "${PLINK}" \
    --bfile kinship_step1 \
    --extract pruned_snps.prune.in \
    --make-bed \
    --out QCed.kinship
fi
append_summary "STEP7_LD_PRUNED" "QCed.kinship" "QCed.kinship"

echo "PCA input QC complete. Output: QCed.kinship.*"
echo "QC summary written to: ${QC_SUMMARY}"

echo "Final outputs:"
echo "  QCed.assoc.*"
echo "  QCed.kinship.*"
echo "  ${QC_SUMMARY}"