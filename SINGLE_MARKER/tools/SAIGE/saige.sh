#!/bin/bash

module load Singularity/4.2.1
module load R/4.2.2

set -euo pipefail

#=============================#
#       INPUT ARGUMENTS      #
#=============================#
if [ $# -lt 5 ]; then
  echo "Usage: $0 <PLINK_DIR> <COVAR_LIST> <QCOVAR_LIST> <MODEL_TAG> <CHR>"
  exit 1
fi

PLINK_DIR=$1
COVAR_LIST=$2
QCOVAR_LIST=$3
MODEL_TAG=$4
CHR=$5

echo "Running SAIGE for chromosome: $CHR"
echo "Covariates: $COVAR_LIST"
echo "Quantitative Covariates: $QCOVAR_LIST"

#=============================#
#         VARIABLES          #
#=============================#
MERGED_PREFIX="${PLINK_DIR}/random_10k_snps"
INPUT_FILE="${PLINK_DIR}/random_10k_snps"
PHENO_FILE="${PLINK_DIR}/pheno.txt"
OUTPUT_DIR="./output"
mkdir -p ${OUTPUT_DIR}/{sparseGRM,saige_output,summary_plots,logs}

#=============================#
#    1. PCA + Pheno with PCs #
#=============================#
PCA_OUT="${OUTPUT_DIR}/mypc"
PHENO_WITH_PCS="${OUTPUT_DIR}/pheno_with_pcs.txt"

../bin/plink --bfile "$MERGED_PREFIX" --pca 10 --out "$PCA_OUT"
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${PCA_OUT}.eigenvec > pcs.tmp
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > pc_title.tmp
cat pc_title.tmp pcs.tmp > pc_ready.tmp
paste "$PHENO_FILE" pc_ready.tmp | sed 's/ /\t/g' > "$PHENO_WITH_PCS"
sed -i '1s/ //g' "$PHENO_WITH_PCS"

#=============================#
#       2. Sparse GRM        #
#=============================#
singularity run Saige_1.3.0.sif createSparseGRM.R \
  --plinkFile="$MERGED_PREFIX" \
  --nThreads=4 \
  --outputPrefix="${OUTPUT_DIR}/sparseGRM/sparseGRM" \
  --numRandomMarkerforSparseKin=2000 \
  --relatednessCutoff=0.125

#=============================#
#     3. Fit Null Model      #
#=============================#
singularity run Saige_1.3.0.sif step1_fitNULLGLMM.R \
  --sparseGRMFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
  --sparseGRMSampleIDFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
  --useSparseGRMtoFitNULL=TRUE \
  --plinkFile="$MERGED_PREFIX" \
  --phenoFile="$PHENO_WITH_PCS" \
  --phenoCol=ADRD \
  --covarColList=${COVAR_LIST} \
  --qCovarColList=${QCOVAR_LIST} \
  --sampleIDColinphenoFile=IID \
  --traitType=binary \
  --outputPrefix="${OUTPUT_DIR}/saige_output/fit_null" \
  --skipVarianceRatioEstimation=FALSE \
  --IsOverwriteVarianceRatioFile=TRUE

#=============================#
#      4. SAIGE Step 2       #
#=============================#
singularity run Saige_1.3.0.sif step2_SPAtests.R \
  --bedFile="${INPUT_FILE}.bed" \
  --bimFile="${INPUT_FILE}.bim" \
  --famFile="${INPUT_FILE}.fam" \
  --AlleleOrder=alt-first \
  --SAIGEOutputFile="${OUTPUT_DIR}/saige_output/saige_results.txt" \
  --GMMATmodelFile="${OUTPUT_DIR}/saige_output/fit_null.rda" \
  --varianceRatioFile="${OUTPUT_DIR}/saige_output/fit_null.varianceRatio.txt" \
  --is_output_moreDetails=TRUE \
  --sparseGRMFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
  --sparseGRMSampleIDFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth=0.05 \
  --LOCO=FALSE \
  --is_fastTest=TRUE

#=============================#
#      5. Post-Processing    #
#=============================#
awk '$7 >= 0.01 && $7 < 0.99' "${OUTPUT_DIR}/saige_output/saige_results.txt" > temp_body

echo "SNP CHR POS Allele1 Allele2 AF BETA SE PVAL N" > title
awk '{print $3,$1,$2,$4,$5,$7,$9,$10,$13,$18+$19}' temp_body > temp1
sed -i '1d' temp1
cat title temp1 | sed 's/ /\t/g' > "${OUTPUT_DIR}/saige_output/sum_stats.txt"

awk '{print $1,$2,$3,$9}' "${OUTPUT_DIR}/saige_output/sum_stats.txt" | sed 's/ /\t/g' > "${OUTPUT_DIR}/saige_output/manhattan_input.txt"

#=============================#
#     6. Manhattan + QQ      #
#=============================#
Rscript - <<EOF
library(CMplot)
library(qqman)

data <- read.table("${OUTPUT_DIR}/saige_output/manhattan_input.txt", header = TRUE, sep = "\t")

# Manhattan plot
png("${OUTPUT_DIR}/summary_plots/manhattan_plot.png", width = 1200, height = 800, res = 150)
CMplot(data,
       plot.type = 'm',
       cex = 1,
       band = 1,
       ylim = c(0, 15),
       col = c("grey30", "grey60"),
       threshold = c(5e-8),
       threshold.col = c("red"),
       threshold.lty = c(5),
       threshold.lwd = c(2),
       amplify = FALSE,
       LOG10 = TRUE)
dev.off()

# QQ plot
p_values <- data\$PVAL
observed_chisq <- qchisq(1 - p_values, df = 1)
lambda <- median(observed_chisq) / 0.4549
lambda_text <- paste0("Lambda = ", round(lambda, 3))

jpeg("${OUTPUT_DIR}/summary_plots/qqplot_with_lambda.jpg", width = 800, height = 800, res = 300)
qq(p_values, main = "QQ Plot of P-values")
text(x = 0.3, y = max(-log10(p_values)) - 0.5, labels = lambda_text, pos = 4, col = "blue", cex = 0.4)
dev.off()
print(paste("Lambda:", lambda))
EOF

#=============================#
#           Cleanup          #
#=============================#
rm -f temp* pcs.tmp pc_title.tmp pc_ready.tmp title
mv *jpg ./output/summary_plots/
mv output/saige_output/sum_stats.txt .
mv ./output/summary_plots/*jpg .

echo "Chromosome ${CHR} processing complete."
