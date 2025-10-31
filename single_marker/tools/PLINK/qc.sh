#!/usr/bin/env bash

input_gwas=$1
input_kinship=$2

# QC for association analysis

# SNP missingness filter
../../../tools/bin/plink \
  --bfile ${input_gwas} \
  --geno 0.02 \
  --make-bed \
  --out assoc_step1

# Individual missingness filter
../../../tools/bin/plink \
  --bfile assoc_step1 \
  --mind 0.02 \
  --make-bed \
  --out assoc_step2

# Minor allele frequency (keep = 0.01, can adjust if you want rare too)
../../../tools/bin/plink \
  --bfile assoc_step2 \
  --maf 0.01 \
  --make-bed \
  --out assoc_step3

# Hardy-Weinberg filter
../../../tools/bin/plink \
  --bfile assoc_step3 \
  --hwe 1e-6 \
  --make-bed \
  --out QCed.assoc

# QC for kinship matrix 

# Start from same input_gwas, remove missing & low MAF SNPs
../../../tools/bin/plink \
  --bfile ${input_kinship} \
  --geno 0.02 \
  --maf 0.05 \
  --mind 0.02 \
  --hwe 1e-6 \
  --make-bed \
  --out kinship_step1

# LD-prune SNPs to avoid correlation
../../../tools/bin/plink \
  --bfile kinship_step1 \
  --indep-pairwise 50 5 0.2 \
  --out pruned_snps

# Keep only pruned SNPs
../../../tools/bin/plink \
  --bfile kinship_step1 \
  --extract pruned_snps.prune.in \
  --make-bed \
  --out QCed.kinship

# cleaning working directory
#rm *_step*
#rm pruned_snps.*

