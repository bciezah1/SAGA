#!/usr/bin/env bash


input_gwas=$1
input_kinship=$2
OUTPUT_DIR=$3
mkdir -p "$OUTPUT_DIR"

./../bin/plink --bfile "$input_gwas" --geno 0.02 --make-bed --out "$OUTPUT_DIR/assoc_step1"
./../bin/plink --bfile "$OUTPUT_DIR/assoc_step1" --mind 0.02 --make-bed --out "$OUTPUT_DIR/assoc_step2"
./../bin/plink --bfile "$OUTPUT_DIR/assoc_step2" --maf 0.01 --make-bed --out "$OUTPUT_DIR/assoc_step3"
./../bin/plink --bfile "$OUTPUT_DIR/assoc_step3" --hwe 1e-6 --make-bed --out "$OUTPUT_DIR/QCed.assoc"

./../bin/plink --bfile "$input_kinship" --geno 0.02 --maf 0.05 --mind 0.02 --hwe 1e-6 --make-bed --out "$OUTPUT_DIR/kinship_step1"
./../bin/plink --bfile "$OUTPUT_DIR/kinship_step1" --indep-pairwise 50 5 0.2 --out "$OUTPUT_DIR/pruned_snps"
./../bin/plink --bfile "$OUTPUT_DIR/kinship_step1" --extract "$OUTPUT_DIR/pruned_snps.prune.in" --make-bed --out "$OUTPUT_DIR/QCed.kinship"
