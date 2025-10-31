#!/usr/bin/env bash
# ===============================================================
# Title: QC PLINK Pipeline for SAIGE and GWAS Preparation
# Author: [Your Name]
# Description:
#   Performs genotype and sample quality control (QC) on PLINK data
#   for both association testing and kinship matrix construction.
#
#   Steps include:
#     1. SNP and sample missingness filters
#     2. MAF and Hardy-Weinberg equilibrium filters
#     3. LD pruning and frequency calculation
#     4. SNP downsampling by MAF bins using an R script
#     5. Preparation of clean datasets for:
#          - Association analysis  ? QCed.assoc.*
#          - Kinship/GRM building  ? QCed.kinship.*
#
# Usage:
#   bash qc_plink_pipeline.sh <input_gwas_prefix> <input_kinship_prefix>
#
# Example:
#   bash qc_plink_pipeline.sh mydata_gwas mydata_kinship
#
# Requirements:
#   - PLINK 1.9 or later
#   - R (with dplyr installed)
#   - R script: step4_downsampling.R
#
# Output:
#   QCed.assoc.*    ? Clean dataset for GWAS/SAIGE Step 1
#   QCed.kinship.*  ? Clean dataset for GRM/sparse GRM
#   list_snps_for_grm.txt ? Downsampled SNPs used for GRM
# ===============================================================

set -e  # Exit immediately if a command exits with a non-zero status



# ----------------------------
# Input arguments
# ----------------------------
input_gwas=$1
input_kinship=$2

if [ -z "$input_gwas" ] || [ -z "$input_kinship" ]; then
  echo "Usage: bash qc_plink_pipeline.sh <input_gwas_prefix> <input_kinship_prefix>"
  exit 1
fi

echo "Starting QC pipeline..."
echo "Input GWAS file prefix: ${input_gwas}"
echo "Input Kinship file prefix: ${input_kinship}"

# ===============================================================
# PART 1: QC FOR ASSOCIATION ANALYSIS
# ===============================================================

echo "Step 1. SNP missingness filter..."
../../bin/plink \
  --bfile ${input_gwas} \
  --geno 0.02 \
  --make-bed \
  --out assoc_step1

echo "Step 2. Individual missingness filter..."
../../bin/plink \
  --bfile assoc_step1 \
  --mind 0.02 \
  --make-bed \
  --out assoc_step2

echo "Step 3. Minor allele frequency filter..."
../../bin/plink \
  --bfile assoc_step2 \
  --maf 0.01 \
  --make-bed \
  --out assoc_step3

echo "Step 4. Hardy-Weinberg equilibrium filter..."
../../bin/plink \
  --bfile assoc_step3 \
  --hwe 1e-6 \
  --make-bed \
  --out QCed.assoc

echo "Association QC complete ? Output: QCed.assoc.*"

# ===============================================================
# PART 2: QC FOR KINSHIP MATRIX / GRM
# ===============================================================

echo "Step 5. LD pruning for kinship..."
../../bin/plink \
  --bfile "${input_kinship}" \
  --indep-pairwise 50 5 0.2 \
  --out temp_snps

echo "Step 6. Extract pruned SNPs..."
../../bin/plink \
  --bfile "${input_kinship}" \
  --extract temp_snps.prune.in \
  --make-bed \
  --out temp_snps_dataset

echo "Step 7. Compute allele frequencies..."
../../bin/plink \
  --bfile temp_snps_dataset \
  --freq \
  --out temp_snps_dataset.freq

echo "Step 8. Downsampling SNPs by MAF using R..."

cp ../step4_downsampling.R .
Rscript step4_downsampling.R

if [ ! -f "list_snps_for_grm.txt" ]; then
  echo "Error: list_snps_for_grm.txt not generated. Check step4_downsampling.R"
  exit 1
fi

echo "Step 9. Extract downsampled SNPs for GRM..."
../../bin/plink \
  --bfile temp_snps_dataset \
  --extract list_snps_for_grm.txt \
  --make-bed \
  --out QCed.kinship

echo "Kinship QC complete ? Output: QCed.kinship.*"

# ===============================================================
# CLEAN-UP
# ===============================================================

echo "Cleaning intermediate files..."
rm -f assoc_step*
rm -f temp_snps.*
rm -f temp_snps_dataset.*
rm -f *.log *.nosex

echo "All steps completed successfully!"
echo "Final outputs:"
echo "  - QCed.assoc.*    ? Association dataset"
echo "  - QCed.kinship.*  ? Kinship dataset"
echo "  - list_snps_for_grm.txt ? SNPs used for GRM"

