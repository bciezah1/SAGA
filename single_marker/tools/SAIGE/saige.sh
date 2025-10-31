#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
module load Singularity/4.2.1
module load R/4.2.2


set -euo pipefail

#=============================#
#       INPUT ARGUMENTS       #
#=============================#
if [ $# -lt 7 ]; then
  echo "Usage: $0 <GENO_INPUT> <PLINK_PHENO> <COVAR_LIST> <QCOVAR_LIST> <PHENO_COL> <TRAIT_TYPE> <WORKING_DIR>"
  echo "TRAIT_TYPE: binary or quantitative"
  exit 1
fi

GENO_INPUT=$1
INPUT_PHENO=$2
COVAR_LIST=$3
QCOVAR_LIST=$4
DIAG_INPUT=$5
TRAIT_TYPE=$6   
WORKING_DIR=$7

#=============================#
#         VARIABLES           #
#=============================#
mkdir -p ${WORKING_DIR}

cd ${WORKING_DIR}
OUTPUT_DIR="./"
mkdir -p ${OUTPUT_DIR}/{sparseGRM,saige_output,summary_plots,logs}

"$SCRIPT_DIR/qc.sh" "../$GENO_INPUT" "../$GENO_INPUT"


echo "Phenotype column: $DIAG_INPUT"
echo "Trait type: $TRAIT_TYPE"
echo "Covariates: $COVAR_LIST"
echo "Quantitative Covariates: $QCOVAR_LIST"

#=============================#
#    1. PCA + Pheno with PCs  #
#=============================#
PCA_OUT="${OUTPUT_DIR}/mypc"
PHENO_WITH_PCS="${OUTPUT_DIR}/pheno_with_pcs.txt"

INPUT_GENO_PCA="QCed.kinship"
INPUT_GENO_KINSHIP="QCed.kinship"
GENO_INPUT="QCed.assoc"

../../bin/plink --bfile "$INPUT_GENO_PCA" --pca 10 --out "$PCA_OUT"
awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${PCA_OUT}.eigenvec > pcs.tmp
echo -e "PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" > pc_title.tmp
cat pc_title.tmp pcs.tmp > pc_ready.tmp
paste "../$INPUT_PHENO" pc_ready.tmp | sed 's/ /\t/g' > "$PHENO_WITH_PCS"
sed -i '1s/ //g' "$PHENO_WITH_PCS"

#=============================#
#       2. Sparse GRM         #
#=============================#
singularity run ../Saige_1.3.0.sif createSparseGRM.R \
  --plinkFile="$INPUT_GENO_KINSHIP" \
  --nThreads=4 \
  --outputPrefix="${OUTPUT_DIR}/sparseGRM/sparseGRM" \
  --numRandomMarkerforSparseKin=2000 \
  --relatednessCutoff=0.125

#=============================#
#     3. Fit Null Model       #
#=============================#
singularity run ../Saige_1.3.0.sif step1_fitNULLGLMM.R \
  --sparseGRMFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx" \
  --sparseGRMSampleIDFile="${OUTPUT_DIR}/sparseGRM/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
  --useSparseGRMtoFitNULL=TRUE \
  --plinkFile="$INPUT_GENO_KINSHIP" \
  --phenoFile="$PHENO_WITH_PCS" \
  --phenoCol=$DIAG_INPUT \
  --covarColList=${COVAR_LIST} \
  --qCovarColList=${QCOVAR_LIST} \
  --sampleIDColinphenoFile=IID \
  --traitType=$TRAIT_TYPE \
  --outputPrefix="${OUTPUT_DIR}/saige_output/fit_null" \
  --skipVarianceRatioEstimation=FALSE \
  --IsOverwriteVarianceRatioFile=TRUE

#=============================#
#      4. SAIGE Step 2        #
#=============================#
singularity run ../Saige_1.3.0.sif step2_SPAtests.R \
  --bedFile="${GENO_INPUT}.bed" \
  --bimFile="${GENO_INPUT}.bim" \
  --famFile="${GENO_INPUT}.fam" \
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
#      5. Post-Processing     #
#=============================#

awk '$7 >= -1 && $7 < 1.2' "${OUTPUT_DIR}/saige_output/saige_results.txt" > temp_body
echo "SNP CHR BP Allele1 Allele2 AF BETA SE PVAL N" > title
awk '{print $3,$1,$2,$4,$5,$7,$9,$10,$13,$18+$19}' temp_body > temp1
sed -i '1d' temp1
cat title temp1 | sed 's/ /\t/g' > "${OUTPUT_DIR}/saige_output/sum_stats.txt"
rm title

awk '{print $1,$2,$3,$9}' "${OUTPUT_DIR}/saige_output/sum_stats.txt" | sed 's/ /\t/g' > "${OUTPUT_DIR}/saige_output/manhattan_input.txt"

#=============================#
#     6. Manhattan + QQ       #
#=============================#

Rscript "$SCRIPT_DIR/create_manhattan.R"  "$OUTPUT_DIR/saige_output/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_qq.plot.R"  "$OUTPUT_DIR/saige_output/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_circular_manhattan.R"  "$OUTPUT_DIR/saige_output/manhattan_input.txt"
Rscript "$SCRIPT_DIR/create_density_plot.R"  "$OUTPUT_DIR/saige_output/manhattan_input.txt"


#=============================#
#           Cleanup           #
#=============================#
# organize
mkdir output 
cd output
mkdir plots tables
cd ..
mv Cir_Manhtn.manhattan_input_circular.jpg circular_manhattan_plot.jpg
mv Rect_Manhtn.manhattan_input_manhattan_highlight.jpg rectangular_manhattan_plot.jpg
mv manhattan_input_qq_lambda.jpg Q_Q_plot.jpg
mv Marker_Density.manhattan_input_density.jpg SNP_density_plot.jpg
mv *jpg ./output/plots
mv ./saige_output/manhattan_input.txt ./output/tables 
mv  pheno_with_pcs.txt ./output/tables
mv ./saige_output/sum_stats.txt ./output/tables
rm -r logs saige_output/ sparseGRM/ summary_plots/ *.tmp tem* list_snps_for_grm.txt step4_downsampling.R
rm mypc.* QCed* 

