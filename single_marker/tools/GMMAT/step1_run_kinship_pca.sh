#!/bin/bash
set -euo pipefail

mkdir -p "$OUTPUT_DIR"

PHENO="$PHENO_INPUT"
PHENO_CLEAN="$OUTPUT_DIR/pheno_for_kinship.txt"
PHENO_GEMMA="$OUTPUT_DIR/pheno_fid_kinship.txt"
PHENO_WITH_PCS="$OUTPUT_DIR/pheno_with_pcs.txt"

echo ">> Running QC on genotype data"
bash qc.sh "$GENO_INPUT" "$GENO_INPUT" "$OUTPUT_DIR"

BFILE="$OUTPUT_DIR/QCed.kinship"
KINSHIP_OUT="mykinship"
PCA_OUT="$OUTPUT_DIR/mypc"

echo ">> Extracting FID, IID, PHENO"
awk '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="FID") fid=i;
    else if ($i=="IID") iid=i;
    else if ($i=="PHENO") adrd=i;
  }
  print "FID IID PHENO"
}
NR>1 {print $fid,$iid,$adrd}
' "$PHENO" > "$PHENO_CLEAN"

tail -n +2 "$PHENO_CLEAN" > "$PHENO_GEMMA"

echo ">> Computing kinship matrix with GEMMA"
echo "THIS IS MY VARIABLE OUTPUT_DIR"
echo "$OUTPUT_DIR"
./../bin/gemma-0.98.5 -bfile "$BFILE" -gk 1 -p "$PHENO_GEMMA" -maf 0.05 -o "$KINSHIP_OUT"
mv ./output/* $OUTPUT_DIR

echo ">> Computing top 10 PCs with PLINK"
./../bin/plink --bfile "$BFILE" --pca 10 --out "$PCA_OUT"

echo ">> Merging PCs with phenotype"
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${PCA_OUT}.eigenvec > "$OUTPUT_DIR/pcs.tmp"
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > "$OUTPUT_DIR/pc_title.tmp"
cat "$OUTPUT_DIR/pc_title.tmp" "$OUTPUT_DIR/pcs.tmp" > "$OUTPUT_DIR/pc_ready.tmp"
echo "oct 22"
echo "$PHENO"
paste "$PHENO" "$OUTPUT_DIR/pc_ready.tmp" | sed 's/ /\t/g' > "$PHENO_WITH_PCS"
