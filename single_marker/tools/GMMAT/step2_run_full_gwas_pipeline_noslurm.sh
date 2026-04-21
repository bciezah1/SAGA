#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <family_type>"
    exit 1
fi

FAMILY_TYPE="$1"
mkdir -p "$OUTPUT_DIR/tmp"
export TMPDIR="$OUTPUT_DIR/tmp"

for chr_number in {1..1}; do
    echo ">> Running GWAS for chromosome $chr_number..."

    Rscript - <<EOF
library(GMMAT)

chr_number <- $chr_number
assoc_prefix <- "$OUTPUT_DIR/QCed.assoc"
model_formula_str <- Sys.getenv("MODEL_FORMULA")

pheno_all <- read.table("$OUTPUT_DIR/pheno_with_pcs.txt", header = TRUE)

if (!("IID" %in% colnames(pheno_all))) {
    stop("pheno_with_pcs.txt must contain an IID column")
}

grm_all <- as.matrix(read.table("$OUTPUT_DIR/mykinship.cXX.txt"))

if (nrow(grm_all) != nrow(pheno_all) || ncol(grm_all) != nrow(pheno_all)) {
    stop("Kinship matrix dimensions do not match phenotype table row count")
}

rownames(grm_all) <- pheno_all\$IID
colnames(grm_all) <- pheno_all\$IID

formula <- as.formula(model_formula_str)
family_type <- tolower("$FAMILY_TYPE")

if (family_type == "quantitative") {
    model0_all <- glmmkin(
        formula,
        kins = grm_all,
        id = "IID",
        data = pheno_all,
        family = gaussian()
    )
} else if (family_type == "binary") {
    model0_all <- glmmkin(
        formula,
        kins = grm_all,
        id = "IID",
        data = pheno_all,
        family = binomial(link = "logit")
    )
} else {
    stop("Invalid family_type. Use 'quantitative' or 'binary'.")
}

outfile_prefix <- paste0("$OUTPUT_DIR/mega_scores_chr", chr_number, "_model1_raw.txt")

glmm.score(
    model0_all,
    infile = assoc_prefix,
    outfile = outfile_prefix,
    nperbatch = 100
)
EOF

    echo ">> Done with chromosome $chr_number"
done