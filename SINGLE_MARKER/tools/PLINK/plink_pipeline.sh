#!/bin/bash

#=============================#
#         SETTINGS           #
#=============================#
set -euo pipefail
module load Plink/1.9.10
module load R/4.2.2

#=============================#
#       INPUT ARGUMENTS      #
#=============================#
if [ $# -lt 3 ]; then
  echo "Usage: $0 <PLINK_DIR> <COVAR_LIST> <QCOVAR_LIST>"
  exit 1
fi

PLINK_DIR=$1
COVAR_LIST=$2
QCOVAR_LIST=$3

MERGED_PREFIX="${PLINK_DIR}/random_10k_snps"
PHENO_FILE="${PLINK_DIR}/pheno.txt"
OUTPUT_DIR="./output"
mkdir -p ${OUTPUT_DIR}/{sparseGRM,saige_output,summary_plots,logs}

echo "Running association analysis"

#=============================#
#    1. PCA + Merge Pheno    #
#=============================#
PCA_OUT="${OUTPUT_DIR}/mypc"
PHENO_WITH_PCS="${OUTPUT_DIR}/pheno_with_pcs.txt"

plink --bfile "$MERGED_PREFIX" --pca 10 --out "$PCA_OUT"

# Add header and merge with phenotype
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${PCA_OUT}.eigenvec > pcs.tmp
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > pc_title.tmp
cat pc_title.tmp pcs.tmp > pc_ready.tmp
paste "$PHENO_FILE" pc_ready.tmp | sed 's/ /\t/g' > "$PHENO_WITH_PCS"
sed -i '1s/ //g' "$PHENO_WITH_PCS"

# Create simplified pheno and covariate files
awk '{print $1,$2,$6}' "$PHENO_WITH_PCS" > pheno
awk '{print $1,$2,$3,$5}' "$PHENO_WITH_PCS" > covariates

#=============================#
#    2. Association Test     #
#=============================#

plink --bfile "$MERGED_PREFIX" \
  --pheno pheno \
  --covar pheno_with_pcs.txt \
  --covar-name SEX,AGE,PC1,PC2,PC3,PC4,PC5 \
  --logistic hide-covar --allow-no-sex --1 --freq \
  --out assoc

#=============================#
#      3. Join Freq + GWAS   #
#=============================#
GWAS_FILE="assoc.assoc.logistic"
FREQ_FILE="assoc.frq"
MERGED_FILE="assoc.merged.txt"

awk '
  NR==FNR && FNR > 1 { freq[$2] = $0; next }
  FNR == 1 {
    print $0 "\tMAF\tNCHROBS"
    next
  }
  {
    split(freq[$2], f, " ")
    print $0 "\t" f[5] "\t" f[6]
  }
' "$FREQ_FILE" "$GWAS_FILE" > "$MERGED_FILE"

sed -i '1d' "$MERGED_FILE"

#=============================#
#     4. Post-Processing     #
#=============================#
awk '$10 >= 0.01 && $10 < 0.99' "$MERGED_FILE" > temp_body
echo -e "CHR SNP BP A1 TEST NMISS OR STAT PVAL MAF NCHROBS" > header.txt
cat header.txt temp_body | sed 's/ /\t/g' > "./sum_stats.txt"

awk '{print $2,$1,$3,$9}' "./sum_stats.txt" | \
  sed 's/ /\t/g' > "./manhattan_input.txt"

#=============================#
#     5. Manhattan + QQ      #
#=============================#
Rscript - <<EOF
library(CMplot)
library(qqman)

data <- read.table("./manhattan_input.txt", header = TRUE, sep = "\t", col.names = c("SNP", "CHR", "BP", "PVAL"))

# Manhattan Plot
png("./manhattan_plot.png", width = 1200, height = 800, res = 150)
CMplot(data,
       plot.type = 'm',
       col = c("grey30", "grey60"),
       LOG10 = TRUE,
       threshold = c(5e-8),
       threshold.col = "red",
       threshold.lty = 5,
       threshold.lwd = 2)
dev.off()

# QQ Plot
p_values <- data\$PVAL
observed_chisq <- qchisq(1 - p_values, df = 1)
lambda <- median(observed_chisq) / 0.4549
lambda_text <- paste0("Lambda = ", round(lambda, 3))

jpeg("./qqplot_with_lambda.jpg", width = 800, height = 800, res = 300)
qq(p_values, main = "QQ Plot of P-values")
text(x = 0.3, y = max(-log10(p_values)) - 0.5, labels = lambda_text, pos = 4, col = "blue", cex = 0.8)
dev.off()
cat("Lambda:", lambda, "\n")
EOF

#=============================#
#         Cleanup (opt)      #
#=============================#
rm -f pcs.tmp pc_title.tmp pc_ready.tmp temp_body header.txt
# Optional: clean logs and intermediate files if needed

echo "processing complete."
