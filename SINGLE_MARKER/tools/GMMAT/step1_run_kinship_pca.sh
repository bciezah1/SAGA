
set -euo pipefail


# File paths
PHENO="${PHENO_INPUT}"
echo "This is my pheno"
echo $PHENO
PHENO_CLEAN=pheno_for_kinship.txt
PHENO_GEMMA=pheno_fid_kinship.txt
BFILE="${KINSHIP_INPUT}"
echo "this is my geno for kinship"
echo $BFILE
echo ">> Using GENO to calculate KINSHIP: $KINSHIP_INPUT"
KINSHIP_OUT=mykinship
PCA_OUT=mypc
PHENO_WITH_PCS=pheno_with_pcs.txt

echo ">> Step 1: Extract FID, IID, ADRD from phenotype file..."
awk '
NR==1 {
  for (i=1; i<=NF; i++) {
    if ($i == "FID") fid=i;
    else if ($i == "IID") iid=i;
    else if ($i == "ADRD") adrd=i;
  }
  print "FID IID ADRD";
}
NR>1 {
  print $fid, $iid, $adrd;
}
' "$PHENO" > "$PHENO_CLEAN"

tail -n +2 "$PHENO_CLEAN" > "$PHENO_GEMMA"

echo ">> Step 2: Compute kinship matrix with GEMMA..."
./../bin/gemma-0.98.5 -bfile "$BFILE" -gk 1 -p "$PHENO_GEMMA" -maf 0.05 -o "$KINSHIP_OUT"

echo ">> Step 3: Compute top 10 principal components with PLINK..."
./../bin/plink --bfile "$BFILE" --pca 10 --out "$PCA_OUT"

echo ">> Step 4: Merge PCs with phenotype file..."
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${PCA_OUT}.eigenvec > pcs.tmp
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > pc_title.tmp
cat pc_title.tmp pcs.tmp > pc_ready.tmp
paste "$PHENO" pc_ready.tmp | sed 's/ /\t/g' > "$PHENO_WITH_PCS"

echo ">> Step 5: Cleaning temporary files..."
rm -f pcs.tmp pc_title.tmp pc_ready.tmp "$PHENO_CLEAN" "$PHENO_GEMMA"
mv ./output/* .

echo ">> Done. Outputs:"
echo "   - Kinship matrix: output/${KINSHIP_OUT}.cXX.txt"
echo "   - PCs: ${PCA_OUT}.eigenvec and ${PCA_OUT}.eigenval"
echo "   - Phenotype with PCs: ${PHENO_WITH_PCS}"
