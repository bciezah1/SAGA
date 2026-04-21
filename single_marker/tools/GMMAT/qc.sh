#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_gwas_prefix> <input_kinship_prefix> <OUTPUT_DIR>"
    exit 1
fi

SCRIPT_DIR="${SAGA_GMMAT_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BIN_DIR="$SCRIPT_DIR/../bin"

input_gwas="$1"
input_kinship="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

"$BIN_DIR/plink" --bfile "$input_gwas" --geno 0.02 --make-bed --out "$OUTPUT_DIR/assoc_step1"
"$BIN_DIR/plink" --bfile "$OUTPUT_DIR/assoc_step1" --mind 0.02 --make-bed --out "$OUTPUT_DIR/assoc_step2"
"$BIN_DIR/plink" --bfile "$OUTPUT_DIR/assoc_step2" --maf 0.01 --make-bed --out "$OUTPUT_DIR/assoc_step3"
"$BIN_DIR/plink" --bfile "$OUTPUT_DIR/assoc_step3" --hwe 1e-6 --make-bed --out "$OUTPUT_DIR/QCed.assoc"

"$BIN_DIR/plink" --bfile "$input_kinship" --geno 0.05 --maf 0.05 --mind 0.15 --hwe 1e-6 --make-bed --out "$OUTPUT_DIR/kinship_step1"
"$BIN_DIR/plink" --bfile "$OUTPUT_DIR/kinship_step1" --indep-pairwise 50 5 0.2 --out "$OUTPUT_DIR/pruned_snps"
"$BIN_DIR/plink" --bfile "$OUTPUT_DIR/kinship_step1" --extract "$OUTPUT_DIR/pruned_snps.prune.in" --make-bed --out "$OUTPUT_DIR/QCed.kinship"