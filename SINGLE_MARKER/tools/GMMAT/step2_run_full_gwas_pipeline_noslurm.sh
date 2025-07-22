#!/bin/bash
set -euo pipefail

# Load R manually if needed (adjust path if R isn't already in PATH)

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
bfile_path <- Sys.getenv("BFILE_PATH_MAF0_05_R2_099")
dosage_dir <- Sys.getenv("DOSAGE_DIR_MAF0_01_R2_080")
model_formula_str <- Sys.getenv("MODEL_FORMULA")

pheno_all <- read.table("./pheno_with_pcs.txt", header=TRUE)
pheno_all\$AGE_normalized <- scale(pheno_all\$AGE)
pheno_all\$id <- seq(1, nrow(pheno_all))

grm_all <- as.matrix(read.table("./mykinship.cXX.txt"))
rownames(grm_all) <- pheno_all\$IID
colnames(grm_all) <- pheno_all\$IID

formula <- as.formula(model_formula_str)
model0_all <- glmmkin(formula, kins = grm_all, id = "IID", data = pheno_all, family = binomial(link = "logit"))

dosage_file <- file.path(dosage_dir, paste0("random_10k_snps"))
prefix <- paste0("mega_scores_chr", chr_number, "_kinship_09_23_2024_5PCs_model1.txt")
glmm.score(model0_all, infile = dosage_file, outfile = prefix, nperbatch = 100)
EOF

  echo ">> Done with chromosome $chr_number"
done
