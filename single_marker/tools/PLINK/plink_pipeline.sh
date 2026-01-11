#!/bin/bash
#===========================================#
#         GWAS PIPELINE (PLINK + R)         #
#===========================================#

set -euo pipefail
SCRIPT_DIR=$(dirname "$(realpath "$0")")
module load R/4.5.2

#===========================================#
#           INPUT ARGUMENTS                 #
#===========================================#
if [ $# -lt 6 ]; then
  echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <COVAR_LIST> <DIAGNOST_VAR> <TRAIT_TYPE> <OUTPUT_DIR>"
  exit 1
fi

GENO_INPUT="$SCRIPT_DIR/$1"       # Input genotype dataset for association
INPUT_PHENO="$SCRIPT_DIR/$2"      # Phenotype file
COVAR_LIST=$3       # Comma-separated list of covariates
DIAGNOST_VAR=$4     # Phenotype/diagnostic variable to test
TRAIT_TYPE=$5       # "binary" or "quantitative"
OUTPUT_DIR=$6       # Output directory for all results

#===========================================#
#           DIRECTORY SETUP                 #
#===========================================#
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "All outputs will be stored in: $OUTPUT_DIR"

#===========================================#
#         TRAIT TYPE SELECTION              #
#===========================================#
if [[ "$TRAIT_TYPE" == "binary" ]]; then
  TEST="--logistic hide-covar --allow-no-sex --1 --freq"
  GWAS_FILE="assoc.assoc.logistic"
elif [[ "$TRAIT_TYPE" == "quantitative" ]]; then
  TEST="--linear hide-covar --allow-no-sex --freq"
  GWAS_FILE="assoc.assoc.linear"
else
  echo "Invalid TRAIT_TYPE: must be 'binary' or 'quantitative'"
  exit 1
fi

echo "Running association analysis for $TRAIT_TYPE trait"

#===========================================#
#        0. QUALITY CONTROL (QC)            #
#===========================================#
# Runs external QC script (produces QCed.assoc / QCed.kinship)
echo "here are my genoe inputs"
echo "$GENO_INPUT"

"$SCRIPT_DIR/qc.sh" "$GENO_INPUT" "$GENO_INPUT" "$OUTPUT_DIR"

GENO_INPUT="./QCed.assoc"
PCA_INPUT="./QCed.kinship"

#===========================================#
#    1. PCA + MERGE WITH PHENOTYPE FILE     #
#===========================================#
PCA_OUT="./mypc"
PHENO_WITH_PCS="./pheno_with_pcs.txt"

# Calculate top 10 PCs
echo "pca input"
echo "$PCA_INPUT"
echo "Current directory: $(pwd)"

../../../tools/bin/plink --bfile "$PCA_INPUT" --pca 10 --out "$PCA_OUT"

# Extract PCs and add header
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "${PCA_OUT}.eigenvec" > "./pcs.tmp"
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > "./pc_title.tmp"
cat "./pc_title.tmp" "./pcs.tmp" > "./pc_ready.tmp"

# Merge PCs with phenotype file
paste "$INPUT_PHENO" "./pc_ready.tmp" | sed 's/ /\t/g' > "$PHENO_WITH_PCS"
sed -i '1s/ //g' "$PHENO_WITH_PCS"

#===========================================#
#         2. ASSOCIATION TESTING            #
#===========================================#
../../../tools/bin/plink --bfile "$GENO_INPUT" \
  --pheno "$PHENO_WITH_PCS" --pheno-name "$DIAGNOST_VAR" \
  --covar "$PHENO_WITH_PCS" --covar-name ${COVAR_LIST} \
  $TEST \
  --out "./assoc"

#===========================================#
#      3. MERGE FREQUENCIES + GWAS          #
#===========================================#
FREQ_FILE="./assoc.frq"
GWAS_FILE="./${GWAS_FILE}"
MERGED_FILE="./assoc.merged.txt"

awk '
  NR==FNR && FNR > 1 { freq[$2] = $0; next }
  FNR == 1 { print $0 "\tMAF\tNCHROBS"; next }
  {
    split(freq[$2], f, " ")
    print $0 "\t" f[5] "\t" f[6]
  }
' "$FREQ_FILE" "$GWAS_FILE" > "$MERGED_FILE"

# Remove duplicated header from merged file
sed -i '1d' "$MERGED_FILE"

#===========================================#
#        4. POST-PROCESSING RESULTS         #
#===========================================#
awk '$10 >= 0.01 && $10 < 0.99' "$MERGED_FILE" > "./temp_body"
echo -e "CHR SNP BP A1 TEST NMISS OR STAT PVAL MAF NCHROBS" > "./header.txt"
cat "./header.txt" "./temp_body" | sed 's/ /\t/g' > "./temporal"

# Clean spacing and save summary stats
awk '{$1=$1; OFS="\t"}1' "./temporal" > "./sum_stats.txt"

# Prepare Manhattan input: SNP, CHR, BP, PVAL
awk '{print $2,$1,$3,$9}' "./sum_stats.txt" | sed 's/ /\t/g' > "./manhattan_input.txt"

#===========================================#
#           5. MANHATTAN + QQ PLOTS         #
#===========================================#

Rscript "$SCRIPT_DIR/create_manhattan.R"  manhattan_input.txt
Rscript "$SCRIPT_DIR/create_qq.plot.R"  manhattan_input.txt
Rscript "$SCRIPT_DIR/create_circular_manhattan.R"  manhattan_input.txt
Rscript "$SCRIPT_DIR/create_density_plot.R"  manhattan_input.txt


#===========================================#
#                CLEANUP                    #
#===========================================#
rm -f "./pcs.tmp" "./pc_title.tmp" "./pc_ready.tmp" \
      "./temp_body" "./header.txt" "./temporal"

echo "Processing complete. Results are in: ${OUTPUT_DIR}"

#===========================================#
#               ORGANIZING                  #
#===========================================#

mkdir output
cd output
mkdir plots tables
cd ../
mv Cir_Manhtn.manhattan_input_circular.jpg circular_manhattan_plot.jpg
mv Rect_Manhtn.manhattan_input_manhattan_highlight.jpg rectangular_manhattan_plot.jpg
mv manhattan_input_qq_lambda.jpg Q_Q_plot.jpg
mv Marker_Density.manhattan_input_density.jpg SNP_density_plot.jpg
mv *jpg ./output/plots
mv manhattan_input.txt  pheno_with_pcs.txt  sum_stats.txt ./output/tables
rm mypc* asso* QCed.* kinship_step1* pruned_snps*


