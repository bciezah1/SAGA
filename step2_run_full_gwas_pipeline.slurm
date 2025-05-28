#!/bin/bash
#SBATCH --job-name=gwas_array
#SBATCH --output=logs/gwas_%A_%a.out
#SBATCH --error=logs/gwas_%A_%a.err
#SBATCH --array=1-22
#SBATCH --time=04:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --exclude=node50,node51,node52,node53,node54,node55,node56,node57,node58,node59,node60,node61,node62,node63,node64,node65,node66,node67,node68,node69

module load R/4.2.2

chr_number=$SLURM_ARRAY_TASK_ID

# Set unique TMPDIR for each task to avoid R_TempDir creation errors
export TMPDIR="$WORKDIR/tmp/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"

cd "$WORKDIR"  # Ensure working directory is set

Rscript - <<EOF
library(GMMAT)

chr_number <- $chr_number
bfile_path <- Sys.getenv("BFILE_PATH_MAF0_05_R2_099")
dosage_dir <- Sys.getenv("DOSAGE_DIR_MAF0_01_R2_080")
model_type <- Sys.getenv("MODEL_TYPE")

pheno_all <- read.table("./pheno_with_pcs.txt", header=TRUE)
pheno_all\$AGE_normalized <- scale(pheno_all\$AGE)
pheno_all\$id <- seq(1, nrow(pheno_all))

grm_all <- as.matrix(read.table("./mykinship.cXX.txt"))
rownames(grm_all) <- pheno_all\$IID
colnames(grm_all) <- pheno_all\$IID

if (model_type == "model1") {
  formula <- as.formula("ADRD ~ AGE + as.factor(SEX) + PC1 + PC2 + PC3 + PC4 + PC5")
} else if (model_type == "model2") {
  formula <- as.formula("ADRD ~ AGE + as.factor(SEX) + PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(APOEe4)")
} else {
  stop("Unknown model type: ", model_type)
}

print(paste("Using formula:", format(formula)))

model0_all <- glmmkin(formula, kins = grm_all, id = "IID", data = pheno_all, family = binomial(link = "logit"))

dosage_file <- file.path(dosage_dir, paste0("chr", chr_number, ".r2.80.maf0.01"))
out_dir <- "./"
prefix <- paste0("mega_scores_chr", chr_number, "_kinship_09_23_2024_5PCs_model1.txt")

glmm.score(model0_all, infile = dosage_file, outfile = paste0(out_dir, prefix), nperbatch = 100)
EOF

# Optional: clean up tmp directory
rm -rf "$TMPDIR"
