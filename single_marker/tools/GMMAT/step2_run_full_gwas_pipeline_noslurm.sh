#!/bin/bash
set -euo pipefail

# ==============================
# Usage:
#   ./run_gmmat.sh <family_type>
#   family_type = "quantitative" ? Continuous (gaussian)
#   family_type = "binary"       ? Binary (binomial logit)
# ==============================

if [ $# -lt 1 ]; then
    echo "Usage: $0 <family_type>"
    echo "  quantitative = Continuous (gaussian)"
    echo "  binary       = Binary (binomial logit)"
    exit 1
fi

FAMILY_TYPE=$1

module load R/4.2.2

cd "$WORKDIR"

for chr_number in {1..1}; do
  echo ">> Running GWAS for chromosome $chr_number..."
  
  TMPDIR="$WORKDIR/tmp/chr_${chr_number}"
  mkdir -p "$TMPDIR"
  export TMPDIR

  Rscript - <<EOF
library(GMMAT)

chr_number <- $chr_number
bfile_path <- Sys.getenv("KINSHIP_INPUT")
dosage_dir <- Sys.getenv("DOSAGE_INPUT")
model_formula_str <- Sys.getenv("MODEL_FORMULA")

cat(">> Kinship file:", bfile_path, "\n")
cat(">> Dosage file:", dosage_dir, "\n")
cat(">> Model formula:", model_formula_str, "\n")

# Load phenotype data
pheno_all <- read.table("./pheno_with_pcs.txt", header=TRUE)
pheno_all\$AGE_normalized <- pheno_all\$AGE
pheno_all\$id <- seq(1, nrow(pheno_all))
print(pheno_all)

# Load kinship matrix
grm_all <- as.matrix(read.table("./mykinship.cXX.txt"))
rownames(grm_all) <- pheno_all\$IID
colnames(grm_all) <- pheno_all\$IID

# Create formula
formula <- as.formula(model_formula_str)

# Select family based on FAMILY_TYPE passed from bash
family_type <- tolower("$FAMILY_TYPE")

if (family_type == "quantitative") {
    cat(">> Running with continuous phenotype (gaussian)\n")
    model0_all <- glmmkin(formula, kins = grm_all, id = "IID", 
                          data = pheno_all, family = gaussian())
} else if (family_type == "binary") {
    cat(">> Running with binary phenotype (binomial logit)\n")
    model0_all <- glmmkin(formula, kins = grm_all, id = "IID", 
                          data = pheno_all, family = binomial(link = "logit"))
} else {
    stop("Invalid family_type: must be 'quantitative' or 'binary'")
}

# Run association
dosage_file <- file.path(dosage_dir)
prefix <- paste0("mega_scores_chr", chr_number, "_kinship_09_23_2024_5PCs_model1.txt")
glmm.score(model0_all, infile = dosage_file, outfile = prefix, nperbatch = 100)
EOF

  echo ">> Done with chromosome $chr_number"
done
