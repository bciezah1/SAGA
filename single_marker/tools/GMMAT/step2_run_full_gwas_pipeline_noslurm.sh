#!/bin/bash
module load R/4.2.2
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <family_type>"
    exit 1
fi

FAMILY_TYPE=$1
mkdir -p "$OUTPUT_DIR/tmp"
export TMPDIR="$OUTPUT_DIR/tmp"

for chr_number in {1..1}; do
    echo ">> Running GWAS for chromosome $chr_number..."

    Rscript - <<EOF
library(GMMAT)

chr_number <- $chr_number
bfile_path <- "$OUTPUT_DIR/QCed.kinship"
dosage_dir <- "$OUTPUT_DIR/QCed.assoc"
model_formula_str <- Sys.getenv("MODEL_FORMULA")

pheno_all <- read.table("$OUTPUT_DIR/pheno_with_pcs.txt", header=TRUE)
cat("First 5 rows and 5 columns:\n")
print(pheno_all[1:5, 1:5])
pheno_all\$AGE_normalized <- pheno_all\$AGE
pheno_all\$id <- seq(1, nrow(pheno_all))

grm_all <- as.matrix(read.table(paste0("$OUTPUT_DIR/mykinship.cXX.txt")))
rownames(grm_all) <- pheno_all\$IID
colnames(grm_all) <- pheno_all\$IID

formula <- as.formula(model_formula_str)
family_type <- tolower("$FAMILY_TYPE")

if(family_type=="quantitative"){
  model0_all <- glmmkin(formula, kins=grm_all, id="IID", data=pheno_all, family=gaussian())
} else if(family_type=="binary"){
  model0_all <- glmmkin(formula, kins=grm_all, id="IID", data=pheno_all, family=binomial(link="logit"))
} else { stop("Invalid family_type") }

prefix <- paste0("$OUTPUT_DIR/mega_scores_chr", chr_number, "_kinship_09_23_2024_5PCs_model1.txt")
glmm.score(model0_all, infile=dosage_dir, outfile=prefix, nperbatch=100)
EOF

    echo ">> Done with chromosome $chr_number"
done
