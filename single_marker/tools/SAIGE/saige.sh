#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if type module >/dev/null 2>&1; then
  module load Singularity/4.2.1 >/dev/null 2>&1 || true
  module load R/4.5.2 >/dev/null 2>&1 || true
fi

PLINK_BIN="${SCRIPT_DIR}/../bin/plink"
QC_SCRIPT="${SCRIPT_DIR}/qc.sh"
SAIGE_IMAGE="${SCRIPT_DIR}/Saige_1.3.0.sif"
MANHATTAN_SCRIPT="${SCRIPT_DIR}/create_manhattan.R"
QQ_SCRIPT="${SCRIPT_DIR}/create_qq.plot.R"
CIRCULAR_SCRIPT="${SCRIPT_DIR}/create_circular_manhattan.R"
DENSITY_SCRIPT="${SCRIPT_DIR}/create_density_plot.R"
SPLIT_SUMSTATS_SCRIPT="${SCRIPT_DIR}/split_sumstats_by_af.R"

ASSOCIATION_MAF_THRESHOLD=${ASSOCIATION_MAF_THRESHOLD:-0.001}
COMMON_RARE_SPLIT_THRESHOLD=${COMMON_RARE_SPLIT_THRESHOLD:-0.01}
CLEANUP_ON_SUCCESS=${CLEANUP_ON_SUCCESS:-true}

usage() {
  cat <<EOF
Usage:
  bash saige.sh <GENO_INPUT> <PHENO_FILE> <COVAR_LIST> <QCOVAR_LIST> <PHENO_COL> <TRAIT_TYPE> <WORKING_DIR>

Arguments:
  GENO_INPUT    PLINK prefix (without .bed/.bim/.fam)
  PHENO_FILE    Phenotype file with header, including FID and IID
  COVAR_LIST    Comma-separated categorical covariates for SAIGE
                Use NONE if not needed
  QCOVAR_LIST   Comma-separated quantitative covariates for SAIGE
                Use NONE if not needed
  PHENO_COL     Phenotype column to analyze
  TRAIT_TYPE    binary or quantitative
  WORKING_DIR   Output working directory

Important SAIGE behavior:
  covarColList must contain ALL covariates in the model
  qCovarColList must be a subset of covarColList containing only quantitative covariates

Examples:
  COVAR_LIST=SEX,COHORT
  QCOVAR_LIST=AGE,PC1,PC2,PC3

  Then SAIGE receives:
    covarColList=SEX,COHORT,AGE,PC1,PC2,PC3
    qCovarColList=AGE,PC1,PC2,PC3

Environment variables:
  ASSOCIATION_MAF_THRESHOLD    Association-branch MAF threshold in qc.sh (default: 0.001)
  COMMON_RARE_SPLIT_THRESHOLD  Split threshold for rare/common output files (currently expected: 0.01)
  CLEANUP_ON_SUCCESS           Remove temporary/intermediate files after successful completion (default: true)
EOF
}

if [ $# -ne 7 ]; then
  usage
  exit 1
fi

GENO_INPUT=$1
INPUT_PHENO=$2
COVAR_LIST_RAW=$3
QCOVAR_LIST_RAW=$4
PHENO_COL=$5
TRAIT_TYPE=$6
WORKING_DIR=$7

if [ ! -x "${PLINK_BIN}" ]; then
  echo "ERROR: PLINK binary not found or not executable: ${PLINK_BIN}"
  exit 1
fi

if [ ! -f "${QC_SCRIPT}" ]; then
  echo "ERROR: QC script not found: ${QC_SCRIPT}"
  exit 1
fi

if [ ! -f "${SAIGE_IMAGE}" ]; then
  echo "ERROR: SAIGE image not found: ${SAIGE_IMAGE}"
  exit 1
fi

if [ ! -f "${SPLIT_SUMSTATS_SCRIPT}" ]; then
  echo "ERROR: Split sumstats script not found: ${SPLIT_SUMSTATS_SCRIPT}"
  exit 1
fi

for f in "${GENO_INPUT}.bed" "${GENO_INPUT}.bim" "${GENO_INPUT}.fam" "${INPUT_PHENO}"; do
  if [ ! -f "${f}" ]; then
    echo "ERROR: Required input file not found: ${f}"
    exit 1
  fi
done

case "${TRAIT_TYPE}" in
  binary|quantitative)
    ;;
  *)
    echo "ERROR: TRAIT_TYPE must be 'binary' or 'quantitative'"
    exit 1
    ;;
esac

if command -v singularity >/dev/null 2>&1; then
  CONTAINER_CMD="singularity"
elif command -v apptainer >/dev/null 2>&1; then
  CONTAINER_CMD="apptainer"
else
  echo "ERROR: Neither singularity nor apptainer was found in PATH."
  exit 1
fi

normalize_csv() {
  local raw="${1:-}"
  if [ -z "${raw}" ]; then
    echo ""
    return 0
  fi
  echo "${raw}" | tr ',' '\n' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | awk '
    NF > 0 &&
    toupper($0) != "NONE" &&
    toupper($0) != "NA" &&
    toupper($0) != "NULL" &&
    !seen[$0]++
  ' | paste -sd, -
}

combine_csv_unique() {
  {
    echo "${1:-}" | tr ',' '\n'
    echo "${2:-}" | tr ',' '\n'
  } | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | awk '
    NF > 0 &&
    toupper($0) != "NONE" &&
    toupper($0) != "NA" &&
    toupper($0) != "NULL" &&
    !seen[$0]++
  ' | paste -sd, -
}

have_file() {
  local f=$1
  [ -s "${f}" ]
}

have_plink_prefix() {
  local prefix=$1
  [ -s "${prefix}.bed" ] && [ -s "${prefix}.bim" ] && [ -s "${prefix}.fam" ]
}

count_data_rows() {
  local f=$1
  if [ ! -f "${f}" ]; then
    echo 0
    return 0
  fi
  awk 'NR>1{n++} END{print n+0}' "${f}"
}

FINAL_CATEGORICAL_COVARS=$(normalize_csv "${COVAR_LIST_RAW}")
FINAL_QUANTITATIVE_COVARS=$(normalize_csv "${QCOVAR_LIST_RAW}")
MODEL_COVARS=$(combine_csv_unique "${FINAL_CATEGORICAL_COVARS}" "${FINAL_QUANTITATIVE_COVARS}")
MODEL_QCOVARS="${FINAL_QUANTITATIVE_COVARS}"

N_THREADS=${N_THREADS:-4}
SPARSE_GRM_RANDOM_MARKERS=${SPARSE_GRM_RANDOM_MARKERS:-2000}
SPARSE_GRM_RELATEDNESS_CUTOFF=${SPARSE_GRM_RELATEDNESS_CUTOFF:-0.125}
PCA_COMPONENTS=${PCA_COMPONENTS:-10}
SPARSE_GRM_MIN_SAMPLES=${SPARSE_GRM_MIN_SAMPLES:-1000}
FORCE_SPARSE_GRM=${FORCE_SPARSE_GRM:-auto}
P_CUTOFF_FOR_FIRTH=${P_CUTOFF_FOR_FIRTH:-0.05}
MAX_CATEGORICAL_LEVEL_RATIO=${MAX_CATEGORICAL_LEVEL_RATIO:-0.20}
MAX_CATEGORICAL_LEVELS_ABS=${MAX_CATEGORICAL_LEVELS_ABS:-20}
MIN_VALID_VARIANCE_RATIO=${MIN_VALID_VARIANCE_RATIO:-1e-12}

mkdir -p "${WORKING_DIR}"
cd "${WORKING_DIR}"

mkdir -p output/plots output/tables logs sparseGRM saige_output intermediate config diagnostics

LOG_FILE="${WORKING_DIR}/logs/pipeline.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=============================================================="
echo "Starting SAIGE pipeline"
echo "Working directory         : ${WORKING_DIR}"
echo "Genotype input            : ${GENO_INPUT}"
echo "Phenotype file            : ${INPUT_PHENO}"
echo "Phenotype column          : ${PHENO_COL}"
echo "Trait type                : ${TRAIT_TYPE}"
echo "Categorical covars input  : ${FINAL_CATEGORICAL_COVARS:-NONE}"
echo "Quantitative covars input : ${FINAL_QUANTITATIVE_COVARS:-NONE}"
echo "SAIGE covarColList        : ${MODEL_COVARS:-NONE}"
echo "SAIGE qCovarColList       : ${MODEL_QCOVARS:-NONE}"
echo "Association MAF threshold : ${ASSOCIATION_MAF_THRESHOLD}"
echo "Rare/common split AF      : ${COMMON_RARE_SPLIT_THRESHOLD}"
echo "Cleanup on success        : ${CLEANUP_ON_SUCCESS}"
echo "Container command         : ${CONTAINER_CMD}"
echo "SAIGE image               : ${SAIGE_IMAGE}"
echo "=============================================================="

{
  echo "RUN_DATE=$(date '+%Y-%m-%d %H:%M:%S')"
  echo "GENO_INPUT=${GENO_INPUT}"
  echo "INPUT_PHENO=${INPUT_PHENO}"
  echo "PHENO_COL=${PHENO_COL}"
  echo "TRAIT_TYPE=${TRAIT_TYPE}"
  echo "FINAL_CATEGORICAL_COVARS=${FINAL_CATEGORICAL_COVARS:-NONE}"
  echo "FINAL_QUANTITATIVE_COVARS=${FINAL_QUANTITATIVE_COVARS:-NONE}"
  echo "MODEL_COVARS=${MODEL_COVARS:-NONE}"
  echo "MODEL_QCOVARS=${MODEL_QCOVARS:-NONE}"
  echo "ASSOCIATION_MAF_THRESHOLD=${ASSOCIATION_MAF_THRESHOLD}"
  echo "COMMON_RARE_SPLIT_THRESHOLD=${COMMON_RARE_SPLIT_THRESHOLD}"
  echo "CLEANUP_ON_SUCCESS=${CLEANUP_ON_SUCCESS}"
  echo "N_THREADS=${N_THREADS}"
  echo "SPARSE_GRM_RANDOM_MARKERS=${SPARSE_GRM_RANDOM_MARKERS}"
  echo "SPARSE_GRM_RELATEDNESS_CUTOFF=${SPARSE_GRM_RELATEDNESS_CUTOFF}"
  echo "PCA_COMPONENTS=${PCA_COMPONENTS}"
  echo "SPARSE_GRM_MIN_SAMPLES=${SPARSE_GRM_MIN_SAMPLES}"
  echo "FORCE_SPARSE_GRM=${FORCE_SPARSE_GRM}"
  echo "P_CUTOFF_FOR_FIRTH=${P_CUTOFF_FOR_FIRTH}"
  echo "MAX_CATEGORICAL_LEVEL_RATIO=${MAX_CATEGORICAL_LEVEL_RATIO}"
  echo "MAX_CATEGORICAL_LEVELS_ABS=${MAX_CATEGORICAL_LEVELS_ABS}"
  echo "MIN_VALID_VARIANCE_RATIO=${MIN_VALID_VARIANCE_RATIO}"
  echo "PLINK_BIN=${PLINK_BIN}"
  echo "QC_SCRIPT=${QC_SCRIPT}"
  echo "SAIGE_IMAGE=${SAIGE_IMAGE}"
  echo "CONTAINER_CMD=${CONTAINER_CMD}"
  echo "R_VERSION=$(R --version | head -n 1)"
  echo "PLINK_VERSION=$("${PLINK_BIN}" --version 2>&1 | head -n 1)"
  echo "CONTAINER_VERSION=$(${CONTAINER_CMD} --version 2>&1 | head -n 1)"
} > config/run_manifest.txt

echo "Step 1. Run QC..."
if have_plink_prefix "QCed.assoc" && have_plink_prefix "QCed.kinship"; then
  echo "  Found existing QCed.assoc.* and QCed.kinship.* ; skipping QC."
else
  bash "${QC_SCRIPT}" "${GENO_INPUT}" "${GENO_INPUT}" "QCed.assoc" "QCed.kinship"
fi

echo "Step 2. Clean phenotype file and standardize missing values..."
CLEANED_PHENO_ALL="intermediate/pheno.cleaned.all.txt"
CLEANED_PHENO_MODEL="intermediate/pheno.analysis_ready.txt"
PHENO_MISSINGNESS_SUMMARY="output/tables/phenotype_missingness_summary.tsv"

if have_file "${CLEANED_PHENO_ALL}" && have_file "${CLEANED_PHENO_MODEL}" && have_file "${PHENO_MISSINGNESS_SUMMARY}"; then
  echo "  Found cleaned phenotype outputs ; skipping."
else
  export SAIGE_INPUT_PHENO="${INPUT_PHENO}"
  export SAIGE_CLEANED_PHENO_ALL="${CLEANED_PHENO_ALL}"
  export SAIGE_CLEANED_PHENO_MODEL="${CLEANED_PHENO_MODEL}"
  export SAIGE_PHENO_MISSINGNESS_SUMMARY="${PHENO_MISSINGNESS_SUMMARY}"
  export SAIGE_PHENO_COL="${PHENO_COL}"
  export SAIGE_TRAIT_TYPE="${TRAIT_TYPE}"
  export SAIGE_MODEL_COVARS="${MODEL_COVARS:-}"
  export SAIGE_MODEL_QCOVARS="${MODEL_QCOVARS:-}"

  Rscript - <<'EOF'
input_pheno <- Sys.getenv("SAIGE_INPUT_PHENO")
output_all <- Sys.getenv("SAIGE_CLEANED_PHENO_ALL")
output_model <- Sys.getenv("SAIGE_CLEANED_PHENO_MODEL")
summary_file <- Sys.getenv("SAIGE_PHENO_MISSINGNESS_SUMMARY")
pheno_col <- Sys.getenv("SAIGE_PHENO_COL")
trait_type <- Sys.getenv("SAIGE_TRAIT_TYPE")
covars_raw <- Sys.getenv("SAIGE_MODEL_COVARS")
qcovars_raw <- Sys.getenv("SAIGE_MODEL_QCOVARS")

csv_to_vec <- function(x) {
  x <- trimws(x)
  if (x == "" || toupper(x) == "NONE") {
    return(character(0))
  }
  vals <- unlist(strsplit(x, ",", fixed = TRUE))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  unique(vals)
}

is_pc_name <- function(x) {
  grepl("^PC([1-9]|10)$", x)
}

covars_all <- csv_to_vec(covars_raw)
qcovars_all <- csv_to_vec(qcovars_raw)

covars_nonpc <- covars_all[!is_pc_name(covars_all)]
qcovars_nonpc <- qcovars_all[!is_pc_name(qcovars_all)]

required_model_cols <- unique(c("FID", "IID", pheno_col, covars_nonpc, qcovars_nonpc))

pheno <- read.table(
  input_pheno,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = "",
  na.strings = c("", "NA", "NaN", "nan", ".", "NULL", "null")
)

if (any(duplicated(colnames(pheno)))) {
  dup_names <- unique(colnames(pheno)[duplicated(colnames(pheno))])
  stop(paste("Duplicate column names found in phenotype file:", paste(dup_names, collapse = ", ")))
}

missing_required <- setdiff(required_model_cols, colnames(pheno))
if (length(missing_required) > 0) {
  stop(paste("Missing required phenotype columns:", paste(missing_required, collapse = ", ")))
}

dup_ids <- duplicated(paste(pheno$FID, pheno$IID, sep = ":"))
if (any(dup_ids)) {
  bad_ids <- unique(paste(pheno$FID[dup_ids], pheno$IID[dup_ids], sep = ":"))
  stop(paste("Duplicate FID/IID pairs found in phenotype file. Example(s):", paste(head(bad_ids, 5), collapse = ", ")))
}

missing_counts <- sapply(required_model_cols, function(col) sum(is.na(pheno[[col]])))
summary_df <- data.frame(
  COLUMN = names(missing_counts),
  N_MISSING = as.integer(missing_counts),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

write.table(
  pheno,
  file = output_all,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

keep_idx <- complete.cases(pheno[, required_model_cols, drop = FALSE])
pheno_model <- pheno[keep_idx, , drop = FALSE]

if (nrow(pheno_model) == 0) {
  stop("No samples remain in phenotype after removing rows with missing required model fields.")
}

if (trait_type == "binary") {
  vals <- sort(unique(pheno_model[[pheno_col]]))
  vals <- vals[!is.na(vals)]
  if (length(vals) != 2) {
    stop(paste("Binary trait must have exactly 2 non-missing values after filtering. Found:", paste(vals, collapse = ", ")))
  }
}

write.table(
  pheno_model,
  file = output_model,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

cat("Phenotype cleaning complete.\n")
cat(paste("All cleaned phenotype written to:", output_all, "\n"))
cat(paste("Analysis-ready phenotype written to:", output_model, "\n"))
cat(paste("Missingness summary written to:", summary_file, "\n"))
EOF
fi

echo "Step 3. Find common participants across association, kinship, and phenotype..."
COMMON_IDS="output/tables/common_participants.txt"
COMMON_COUNTS="output/tables/common_participant_counts.tsv"

if have_file "${COMMON_IDS}" && have_file "${COMMON_COUNTS}"; then
  echo "  Found common participant outputs ; skipping."
else
  export SAIGE_COMMON_IDS="${COMMON_IDS}"
  export SAIGE_COMMON_COUNTS="${COMMON_COUNTS}"
  export SAIGE_CLEANED_PHENO_MODEL_PATH="${CLEANED_PHENO_MODEL}"

  Rscript - <<'EOF'
assoc_fam <- read.table("QCed.assoc.fam", header = FALSE, stringsAsFactors = FALSE)
kin_fam <- read.table("QCed.kinship.fam", header = FALSE, stringsAsFactors = FALSE)
pheno <- read.table(Sys.getenv("SAIGE_CLEANED_PHENO_MODEL_PATH"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

colnames(assoc_fam)[1:2] <- c("FID", "IID")
colnames(kin_fam)[1:2] <- c("FID", "IID")

assoc_fam <- assoc_fam[, c("FID", "IID")]
kin_fam <- kin_fam[, c("FID", "IID")]
pheno <- pheno[, c("FID", "IID")]

assoc_fam$key <- paste(assoc_fam$FID, assoc_fam$IID, sep = ":")
kin_fam$key <- paste(kin_fam$FID, kin_fam$IID, sep = ":")
pheno$key <- paste(pheno$FID, pheno$IID, sep = ":")

common_keys <- Reduce(intersect, list(assoc_fam$key, kin_fam$key, pheno$key))

if (length(common_keys) == 0) {
  stop("No common participants found across QCed.assoc, QCed.kinship, and phenotype.")
}

common_kin_order <- kin_fam[kin_fam$key %in% common_keys, c("FID", "IID")]

write.table(
  common_kin_order,
  file = Sys.getenv("SAIGE_COMMON_IDS"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

counts <- data.frame(
  SOURCE = c("QCed.assoc", "QCed.kinship", "pheno.analysis_ready", "common_participants"),
  N = c(nrow(assoc_fam), nrow(kin_fam), nrow(pheno), nrow(common_kin_order)),
  stringsAsFactors = FALSE
)

write.table(
  counts,
  file = Sys.getenv("SAIGE_COMMON_COUNTS"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Common participant list created successfully.\n")
EOF
fi

echo "Step 4. Restrict association and kinship datasets to common participants..."
if have_plink_prefix "matched.assoc" && have_plink_prefix "matched.kinship"; then
  echo "  Found matched.assoc.* and matched.kinship.* ; skipping."
else
  "${PLINK_BIN}" \
    --bfile QCed.assoc \
    --keep "${COMMON_IDS}" \
    --make-bed \
    --out matched.assoc

  "${PLINK_BIN}" \
    --bfile QCed.kinship \
    --keep "${COMMON_IDS}" \
    --make-bed \
    --out matched.kinship
fi

N_ASSOC=$(wc -l < matched.assoc.fam)
N_KINSHIP=$(wc -l < matched.kinship.fam)
N_COMMON=$(wc -l < "${COMMON_IDS}")

if [ "${N_ASSOC}" -ne "${N_COMMON}" ] || [ "${N_KINSHIP}" -ne "${N_COMMON}" ]; then
  echo "ERROR: Sample count mismatch after restricting to common participants."
  echo "  common_ids      : ${N_COMMON}"
  echo "  matched.assoc   : ${N_ASSOC}"
  echo "  matched.kinship : ${N_KINSHIP}"
  exit 1
fi

echo "Matched sample count: ${N_COMMON}"

if [ "${FORCE_SPARSE_GRM}" = "true" ]; then
  USE_SPARSE_GRM="true"
elif [ "${FORCE_SPARSE_GRM}" = "false" ]; then
  USE_SPARSE_GRM="false"
else
  if [ "${N_COMMON}" -ge "${SPARSE_GRM_MIN_SAMPLES}" ]; then
    USE_SPARSE_GRM="true"
  else
    USE_SPARSE_GRM="false"
  fi
fi

echo "Sparse GRM decision:"
echo "  FORCE_SPARSE_GRM       : ${FORCE_SPARSE_GRM}"
echo "  SPARSE_GRM_MIN_SAMPLES : ${SPARSE_GRM_MIN_SAMPLES}"
echo "  matched samples        : ${N_COMMON}"
echo "  USE_SPARSE_GRM         : ${USE_SPARSE_GRM}"

echo "Step 5. Run PCA on matched kinship dataset..."
PCA_OUT="intermediate/mypc"
if have_file "${PCA_OUT}.eigenvec" && have_file "${PCA_OUT}.eigenval"; then
  echo "  Found PCA outputs ; skipping."
else
  "${PLINK_BIN}" \
    --bfile matched.kinship \
    --pca "${PCA_COMPONENTS}" \
    --out "${PCA_OUT}"
fi

echo "Step 6. Merge phenotype with PCs and align to matched.kinship sample order..."
PHENO_WITH_PCS="output/tables/pheno_with_pcs.txt"

if have_file "${PHENO_WITH_PCS}"; then
  echo "  Found pheno_with_pcs.txt ; skipping."
else
  export SAIGE_PHENO_WITH_PCS="${PHENO_WITH_PCS}"
  export SAIGE_PCA_OUT="${PCA_OUT}"
  export SAIGE_PCA_COMPONENTS="${PCA_COMPONENTS}"
  export SAIGE_CLEANED_PHENO_MODEL_PATH="${CLEANED_PHENO_MODEL}"

  Rscript - <<'EOF'
pheno <- read.table(Sys.getenv("SAIGE_CLEANED_PHENO_MODEL_PATH"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
pca <- read.table(paste0(Sys.getenv("SAIGE_PCA_OUT"), ".eigenvec"), header = FALSE, stringsAsFactors = FALSE)
fam <- read.table("matched.kinship.fam", header = FALSE, stringsAsFactors = FALSE)

n_pcs <- as.integer(Sys.getenv("SAIGE_PCA_COMPONENTS"))
expected_cols <- 2 + n_pcs

if (ncol(pca) < expected_cols) {
  stop(paste("Unexpected PCA eigenvec format. Expected at least", expected_cols, "columns."))
}

colnames(pca) <- c("FID", "IID", paste0("PC", seq_len(n_pcs)))
fam <- fam[, 1:2]
colnames(fam) <- c("FID", "IID")

merged <- merge(pheno, pca, by = c("FID", "IID"), all = FALSE)
merged$key <- paste(merged$FID, merged$IID, sep = ":")
fam$key <- paste(fam$FID, fam$IID, sep = ":")

missing_pcs <- setdiff(fam$key, merged$key)
if (length(missing_pcs) > 0) {
  stop("Some matched kinship participants are missing after phenotype-PCA merge.")
}

merged <- merged[match(fam$key, merged$key), ]
merged$key <- NULL

write.table(
  merged,
  file = Sys.getenv("SAIGE_PHENO_WITH_PCS"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

cat("Phenotype with PCs created successfully.\n")
EOF
fi

echo "Step 7. Covariate sanity checks on analysis-ready phenotype with PCs..."
COVARIATE_SANITY_REPORT="diagnostics/covariate_sanity_report.tsv"

if have_file "${COVARIATE_SANITY_REPORT}"; then
  echo "  Found covariate sanity report ; skipping."
else
  export SAIGE_COVARIATE_SANITY_REPORT="${COVARIATE_SANITY_REPORT}"
  export SAIGE_MAX_CATEGORICAL_LEVEL_RATIO="${MAX_CATEGORICAL_LEVEL_RATIO}"
  export SAIGE_MAX_CATEGORICAL_LEVELS_ABS="${MAX_CATEGORICAL_LEVELS_ABS}"
  export SAIGE_MODEL_COVARS="${MODEL_COVARS:-}"
  export SAIGE_MODEL_QCOVARS="${MODEL_QCOVARS:-}"
  export SAIGE_PHENO_WITH_PCS="${PHENO_WITH_PCS}"

  Rscript - <<'EOF'
pheno <- read.table(Sys.getenv("SAIGE_PHENO_WITH_PCS"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

csv_to_vec <- function(x) {
  x <- trimws(x)
  if (x == "" || toupper(x) == "NONE") {
    return(character(0))
  }
  vals <- unlist(strsplit(x, ",", fixed = TRUE))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  unique(vals)
}

model_covars <- csv_to_vec(Sys.getenv("SAIGE_MODEL_COVARS"))
model_qcovars <- csv_to_vec(Sys.getenv("SAIGE_MODEL_QCOVARS"))
categorical_covars <- setdiff(model_covars, model_qcovars)

max_level_ratio <- as.numeric(Sys.getenv("SAIGE_MAX_CATEGORICAL_LEVEL_RATIO"))
max_levels_abs <- as.integer(Sys.getenv("SAIGE_MAX_CATEGORICAL_LEVELS_ABS"))
n <- nrow(pheno)

report <- data.frame(
  VARIABLE = character(0),
  TYPE = character(0),
  N_UNIQUE = integer(0),
  STATUS = character(0),
  DETAIL = character(0),
  stringsAsFactors = FALSE
)

append_row <- function(var, type, n_unique, status, detail) {
  data.frame(
    VARIABLE = var,
    TYPE = type,
    N_UNIQUE = n_unique,
    STATUS = status,
    DETAIL = detail,
    stringsAsFactors = FALSE
  )
}

for (v in categorical_covars) {
  if (!v %in% colnames(pheno)) {
    stop(paste("Categorical covariate missing from pheno_with_pcs:", v))
  }
  vals <- pheno[[v]]
  vals <- vals[!is.na(vals)]
  n_unique <- length(unique(vals))

  if (n_unique <= 1) {
    report <- rbind(report, append_row(v, "categorical", n_unique, "FAIL", "Covariate is constant after filtering"))
    next
  }

  if (n_unique > max_levels_abs || n_unique > ceiling(n * max_level_ratio)) {
    report <- rbind(report, append_row(v, "categorical", n_unique, "WARN", "Too many levels for cohort size; may destabilize the model"))
  } else {
    report <- rbind(report, append_row(v, "categorical", n_unique, "OK", ""))
  }
}

for (v in model_qcovars) {
  if (!v %in% colnames(pheno)) {
    stop(paste("Quantitative covariate missing from pheno_with_pcs:", v))
  }
  vals <- pheno[[v]]
  vals <- vals[!is.na(vals)]
  n_unique <- length(unique(vals))

  if (n_unique <= 1) {
    report <- rbind(report, append_row(v, "quantitative", n_unique, "FAIL", "Covariate is constant after filtering"))
  } else {
    report <- rbind(report, append_row(v, "quantitative", n_unique, "OK", ""))
  }
}

write.table(
  report,
  file = Sys.getenv("SAIGE_COVARIATE_SANITY_REPORT"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

if (any(report$STATUS == "FAIL")) {
  bad <- report$VARIABLE[report$STATUS == "FAIL"]
  stop(paste("Covariate sanity checks failed for:", paste(bad, collapse = ", ")))
}

cat("Covariate sanity checks completed.\n")
EOF
fi

SPARSE_GRM_FILE=""
SPARSE_GRM_IDS=""

if [ "${USE_SPARSE_GRM}" = "true" ]; then
  echo "Step 8. Build sparse GRM..."
  SPARSE_GRM_FILE="sparseGRM/sparseGRM_relatednessCutoff_${SPARSE_GRM_RELATEDNESS_CUTOFF}_${SPARSE_GRM_RANDOM_MARKERS}_randomMarkersUsed.sparseGRM.mtx"
  SPARSE_GRM_IDS="sparseGRM/sparseGRM_relatednessCutoff_${SPARSE_GRM_RELATEDNESS_CUTOFF}_${SPARSE_GRM_RANDOM_MARKERS}_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

  if have_file "${SPARSE_GRM_FILE}" && have_file "${SPARSE_GRM_IDS}"; then
    echo "  Found sparse GRM outputs ; skipping."
  else
    "${CONTAINER_CMD}" run "${SAIGE_IMAGE}" createSparseGRM.R \
      --plinkFile="matched.kinship" \
      --nThreads="${N_THREADS}" \
      --outputPrefix="sparseGRM/sparseGRM" \
      --numRandomMarkerforSparseKin="${SPARSE_GRM_RANDOM_MARKERS}" \
      --relatednessCutoff="${SPARSE_GRM_RELATEDNESS_CUTOFF}"
  fi

  if [ ! -f "${SPARSE_GRM_FILE}" ] || [ ! -f "${SPARSE_GRM_IDS}" ]; then
    echo "ERROR: Sparse GRM output files were not created as expected."
    exit 1
  fi
else
  echo "Step 8. Sparse GRM skipped for this cohort size."
fi

echo "Step 9. Fit SAIGE null model..."
FIT_NULL_RDA="saige_output/fit_null.rda"
FIT_NULL_VR="saige_output/fit_null.varianceRatio.txt"

if have_file "${FIT_NULL_RDA}" && have_file "${FIT_NULL_VR}"; then
  echo "  Found null model outputs ; skipping."
else
  STEP1_ARGS=(
    --plinkFile="matched.kinship"
    --phenoFile="${PHENO_WITH_PCS}"
    --phenoCol="${PHENO_COL}"
    --sampleIDColinphenoFile=IID
    --traitType="${TRAIT_TYPE}"
    --outputPrefix="saige_output/fit_null"
    --skipVarianceRatioEstimation=FALSE
    --IsOverwriteVarianceRatioFile=TRUE
    --nThreads="${N_THREADS}"
  )

  if [ -n "${MODEL_COVARS:-}" ]; then
    STEP1_ARGS+=(--covarColList="${MODEL_COVARS}")
  fi

  if [ -n "${MODEL_QCOVARS:-}" ]; then
    STEP1_ARGS+=(--qCovarColList="${MODEL_QCOVARS}")
  fi

  if [ "${USE_SPARSE_GRM}" = "true" ]; then
    STEP1_ARGS+=(
      --useSparseGRMtoFitNULL=TRUE
      --sparseGRMFile="${SPARSE_GRM_FILE}"
      --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}"
    )
  else
    STEP1_ARGS+=(--useSparseGRMtoFitNULL=FALSE)
  fi

  "${CONTAINER_CMD}" run "${SAIGE_IMAGE}" step1_fitNULLGLMM.R "${STEP1_ARGS[@]}"
fi

if [ ! -f "${FIT_NULL_RDA}" ] || [ ! -f "${FIT_NULL_VR}" ]; then
  echo "ERROR: Null model outputs were not created."
  exit 1
fi

echo "Step 10. Validate null-model outputs..."
NULL_MODEL_DIAGNOSTICS="diagnostics/null_model_diagnostics.tsv"

if have_file "${NULL_MODEL_DIAGNOSTICS}"; then
  echo "  Found null model diagnostics ; skipping."
else
  export SAIGE_NULL_MODEL_DIAGNOSTICS="${NULL_MODEL_DIAGNOSTICS}"
  export SAIGE_MIN_VALID_VARIANCE_RATIO="${MIN_VALID_VARIANCE_RATIO}"

  Rscript - <<'EOF'
vr_file <- "saige_output/fit_null.varianceRatio.txt"
out_file <- Sys.getenv("SAIGE_NULL_MODEL_DIAGNOSTICS")
min_valid_vr <- as.numeric(Sys.getenv("SAIGE_MIN_VALID_VARIANCE_RATIO"))

vr <- read.table(vr_file, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
if (nrow(vr) == 0) {
  stop("Variance ratio file is empty.")
}

first_col <- suppressWarnings(as.numeric(vr[[1]]))
if (all(is.na(first_col))) {
  stop("Could not parse variance ratio values from fit_null.varianceRatio.txt")
}

status <- rep("OK", length(first_col))
detail <- rep("", length(first_col))

for (i in seq_along(first_col)) {
  val <- first_col[i]
  if (is.na(val) || !is.finite(val)) {
    status[i] <- "FAIL"
    detail[i] <- "Variance ratio is NA or non-finite"
  } else if (val <= 0) {
    status[i] <- "FAIL"
    detail[i] <- "Variance ratio is non-positive"
  } else if (val < min_valid_vr) {
    status[i] <- "FAIL"
    detail[i] <- paste("Variance ratio is too close to zero; threshold =", min_valid_vr)
  }
}

diagnostics <- data.frame(
  VARIANCE_RATIO = first_col,
  LABEL = if (ncol(vr) >= 2) vr[[2]] else "",
  EXTRA = if (ncol(vr) >= 3) vr[[3]] else "",
  STATUS = status,
  DETAIL = detail,
  stringsAsFactors = FALSE
)

write.table(
  diagnostics,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

if (any(diagnostics$STATUS == "FAIL")) {
  print(diagnostics)
  stop("Null model validation failed. See diagnostics/null_model_diagnostics.tsv")
}

cat("Null model validation passed.\n")
EOF
fi

echo "Variance ratio file contents:"
cat "saige_output/fit_null.varianceRatio.txt"

echo "Step 11. Run SAIGE step 2 association tests..."
SAIGE_RESULTS="saige_output/saige_results.txt"

if have_file "${SAIGE_RESULTS}"; then
  echo "  Found saige_results.txt ; skipping."
else
  STEP2_ARGS=(
    --bedFile="matched.assoc.bed"
    --bimFile="matched.assoc.bim"
    --famFile="matched.assoc.fam"
    --AlleleOrder=alt-first
    --SAIGEOutputFile="${SAIGE_RESULTS}"
    --GMMATmodelFile="${FIT_NULL_RDA}"
    --varianceRatioFile="${FIT_NULL_VR}"
    --is_output_moreDetails=TRUE
    --is_Firth_beta=TRUE
    --pCutoffforFirth="${P_CUTOFF_FOR_FIRTH}"
    --LOCO=FALSE
    --is_fastTest=TRUE
  )

  if [ "${USE_SPARSE_GRM}" = "true" ]; then
    STEP2_ARGS+=(
      --sparseGRMFile="${SPARSE_GRM_FILE}"
      --sparseGRMSampleIDFile="${SPARSE_GRM_IDS}"
    )
  fi

  "${CONTAINER_CMD}" run "${SAIGE_IMAGE}" step2_SPAtests.R "${STEP2_ARGS[@]}"
fi

if [ ! -f "${SAIGE_RESULTS}" ]; then
  echo "ERROR: SAIGE results file was not created."
  exit 1
fi

echo "Step 12. Post-process SAIGE results into summary statistics..."
SUM_STATS="output/tables/sum_stats.txt"
MANHATTAN_INPUT="output/tables/manhattan_input.txt"

if have_file "${SUM_STATS}" && have_file "${MANHATTAN_INPUT}"; then
  echo "  Found sum_stats.txt and manhattan_input.txt ; skipping."
else
  export SAIGE_SUM_STATS="${SUM_STATS}"
  export SAIGE_MANHATTAN_INPUT="${MANHATTAN_INPUT}"

  Rscript - <<'EOF'
pick_col <- function(df, candidates, required = TRUE) {
  idx <- which(colnames(df) %in% candidates)
  if (length(idx) == 0) {
    if (required) {
      stop(paste("Could not find any of the required columns:", paste(candidates, collapse = ", ")))
    } else {
      return(NULL)
    }
  }
  colnames(df)[idx[1]]
}

res <- read.table(
  "saige_output/saige_results.txt",
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

chr_col <- pick_col(res, c("CHR", "Chr", "chr"))
bp_col <- pick_col(res, c("POS", "BP", "Pos", "bp"))
snp_col <- pick_col(res, c("MarkerID", "SNPID", "SNP", "markerID", "ID"))
a1_col <- pick_col(res, c("Allele1", "A1", "Allele_1"))
a2_col <- pick_col(res, c("Allele2", "A2", "Allele_2"))
af_col <- pick_col(res, c("AF_Allele2", "AF", "AltAF", "AlleleFreq"), required = FALSE)
beta_col <- pick_col(res, c("BETA", "Beta", "beta"))
se_col <- pick_col(res, c("SE", "se"))
p_col <- pick_col(res, c("p.value", "PVAL", "pval", "P", "p_value"))

n_col <- pick_col(res, c("N", "n"), required = FALSE)
n_case_col <- pick_col(res, c("N_case", "N.Cases", "NCase"), required = FALSE)
n_ctrl_col <- pick_col(res, c("N_ctrl", "N.Controls", "NCtrl"), required = FALSE)

if (!is.null(n_col)) {
  n_vec <- res[[n_col]]
} else if (!is.null(n_case_col) && !is.null(n_ctrl_col)) {
  n_vec <- res[[n_case_col]] + res[[n_ctrl_col]]
} else {
  n_vec <- NA
}

if (!is.null(af_col)) {
  af_vec <- res[[af_col]]
} else {
  af_vec <- NA
}

sumstats <- data.frame(
  SNP = res[[snp_col]],
  CHR = res[[chr_col]],
  BP = res[[bp_col]],
  Allele1 = res[[a1_col]],
  Allele2 = res[[a2_col]],
  AF = af_vec,
  BETA = res[[beta_col]],
  SE = res[[se_col]],
  PVAL = res[[p_col]],
  N = n_vec,
  stringsAsFactors = FALSE
)

sumstats <- sumstats[!is.na(sumstats$SNP) & !is.na(sumstats$CHR) & !is.na(sumstats$BP) & !is.na(sumstats$PVAL), ]

if (!all(is.na(sumstats$AF))) {
  sumstats <- sumstats[is.na(sumstats$AF) | (sumstats$AF >= 0 & sumstats$AF <= 1), ]
}

write.table(
  sumstats,
  file = Sys.getenv("SAIGE_SUM_STATS"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

manhattan_df <- sumstats[, c("CHR", "BP", "SNP", "PVAL")]
write.table(
  manhattan_df,
  file = Sys.getenv("SAIGE_MANHATTAN_INPUT"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

cat("Summary statistics and Manhattan input written successfully.\n")
EOF
fi

echo "Step 12b. Split summary statistics into rare and common outputs..."
SUM_STATS_RARE="output/tables/sum_stats_AF_lt_0.01.txt"
SUM_STATS_COMMON="output/tables/sum_stats_AF_ge_0.01.txt"
MANHATTAN_INPUT_RARE="output/tables/sum_stats_AF_lt_0.01_CHR_BP_SNP_PVAL.txt"
MANHATTAN_INPUT_COMMON="output/tables/sum_stats_AF_ge_0.01_CHR_BP_SNP_PVAL.txt"

if have_file "${SUM_STATS_RARE}" && have_file "${SUM_STATS_COMMON}" && have_file "${MANHATTAN_INPUT_RARE}" && have_file "${MANHATTAN_INPUT_COMMON}"; then
  echo "  Found split rare/common outputs ; skipping."
else
  if [ "${COMMON_RARE_SPLIT_THRESHOLD}" != "0.01" ]; then
    echo "WARNING: split_sumstats_by_af.R currently writes files using the built-in 0.01 threshold naming."
    echo "WARNING: Continuing with the existing rare/common split script."
  fi

  Rscript "${SPLIT_SUMSTATS_SCRIPT}" "${SUM_STATS}"
fi

echo "Step 13. Create rare/common plot sets directly in output/plots..."

COMMON_PLOT_INPUT="output/plots/common_plot_input.txt"
RARE_PLOT_INPUT="output/plots/rare_plot_input.txt"

COMMON_MANHATTAN_PLOT="output/plots/Rect_Manhtn.common_plot_input.jpg"
RARE_MANHATTAN_PLOT="output/plots/Rect_Manhtn.rare_plot_input.jpg"

COMMON_QQ_PLOT="output/plots/common_plot_input_qq_lambda.jpg"
RARE_QQ_PLOT="output/plots/rare_plot_input_qq_lambda.jpg"

COMMON_CIRCULAR_PLOT="output/plots/common_plot_input_circular.jpg"
RARE_CIRCULAR_PLOT="output/plots/rare_plot_input_circular.jpg"

COMMON_DENSITY_PLOT="output/plots/common_plot_input_density.jpg"
RARE_DENSITY_PLOT="output/plots/rare_plot_input_density.jpg"

COMMON_SKIP_NOTE="output/plots/common_plots_skipped.txt"
RARE_SKIP_NOTE="output/plots/rare_plots_skipped.txt"

make_skip_note() {
  local note_file=$1
  local label=$2
  local reason=$3

  {
    echo "${label} plot set skipped"
    echo "Reason: ${reason}"
    echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
  } > "${note_file}"
}

prepare_plot_input() {
  local source_input=$1
  local dest_input=$2
  local label=$3
  local skip_note=$4
  local n_rows

  n_rows=$(count_data_rows "${source_input}")
  if [ "${n_rows}" -eq 0 ]; then
    echo "  ${label}: no variants found; skipping plot set."
    make_skip_note "${skip_note}" "${label}" "No variants available in ${source_input}"
    return 1
  fi

  cp -f "${source_input}" "${dest_input}"
  return 0
}

make_common_plots=false
make_rare_plots=false

if prepare_plot_input "${MANHATTAN_INPUT_COMMON}" "${COMMON_PLOT_INPUT}" "common" "${COMMON_SKIP_NOTE}"; then
  make_common_plots=true
fi

if prepare_plot_input "${MANHATTAN_INPUT_RARE}" "${RARE_PLOT_INPUT}" "rare" "${RARE_SKIP_NOTE}"; then
  make_rare_plots=true
fi

if [ "${make_common_plots}" = true ]; then
  if [ ! -f "${MANHATTAN_SCRIPT}" ]; then
    echo "ERROR: Manhattan script not found: ${MANHATTAN_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${QQ_SCRIPT}" ]; then
    echo "ERROR: QQ script not found: ${QQ_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${CIRCULAR_SCRIPT}" ]; then
    echo "ERROR: Circular Manhattan script not found: ${CIRCULAR_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${DENSITY_SCRIPT}" ]; then
    echo "ERROR: Density script not found: ${DENSITY_SCRIPT}"
    exit 1
  fi

  if have_file "${COMMON_MANHATTAN_PLOT}"; then
    echo "  Found $(basename "${COMMON_MANHATTAN_PLOT}") ; skipping."
  else
    Rscript "${MANHATTAN_SCRIPT}" "${COMMON_PLOT_INPUT}"
  fi

  if have_file "${COMMON_QQ_PLOT}"; then
    echo "  Found $(basename "${COMMON_QQ_PLOT}") ; skipping."
  else
    Rscript "${QQ_SCRIPT}" "${COMMON_PLOT_INPUT}" "output/plots"
  fi

  if have_file "${COMMON_CIRCULAR_PLOT}"; then
    echo "  Found $(basename "${COMMON_CIRCULAR_PLOT}") ; skipping."
  else
    Rscript "${CIRCULAR_SCRIPT}" "${COMMON_PLOT_INPUT}" "output/plots"
  fi

  if have_file "${COMMON_DENSITY_PLOT}"; then
    echo "  Found $(basename "${COMMON_DENSITY_PLOT}") ; skipping."
  else
    Rscript "${DENSITY_SCRIPT}" "${COMMON_PLOT_INPUT}" "output/plots"
  fi
fi

if [ "${make_rare_plots}" = true ]; then
  if [ ! -f "${MANHATTAN_SCRIPT}" ]; then
    echo "ERROR: Manhattan script not found: ${MANHATTAN_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${QQ_SCRIPT}" ]; then
    echo "ERROR: QQ script not found: ${QQ_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${CIRCULAR_SCRIPT}" ]; then
    echo "ERROR: Circular Manhattan script not found: ${CIRCULAR_SCRIPT}"
    exit 1
  fi
  if [ ! -f "${DENSITY_SCRIPT}" ]; then
    echo "ERROR: Density script not found: ${DENSITY_SCRIPT}"
    exit 1
  fi

  if have_file "${RARE_MANHATTAN_PLOT}"; then
    echo "  Found $(basename "${RARE_MANHATTAN_PLOT}") ; skipping."
  else
    Rscript "${MANHATTAN_SCRIPT}" "${RARE_PLOT_INPUT}"
  fi

  if have_file "${RARE_QQ_PLOT}"; then
    echo "  Found $(basename "${RARE_QQ_PLOT}") ; skipping."
  else
    Rscript "${QQ_SCRIPT}" "${RARE_PLOT_INPUT}" "output/plots"
  fi

  if have_file "${RARE_CIRCULAR_PLOT}"; then
    echo "  Found $(basename "${RARE_CIRCULAR_PLOT}") ; skipping."
  else
    Rscript "${CIRCULAR_SCRIPT}" "${RARE_PLOT_INPUT}" "output/plots"
  fi

  if have_file "${RARE_DENSITY_PLOT}"; then
    echo "  Found $(basename "${RARE_DENSITY_PLOT}") ; skipping."
  else
    Rscript "${DENSITY_SCRIPT}" "${RARE_PLOT_INPUT}" "output/plots"
  fi
fi

echo "Step 14. Write final run summary..."
{
  echo "SAIGE pipeline completed successfully."
  echo
  echo "Key outputs:"
  echo "  output/tables/common_participants.txt"
  echo "  output/tables/common_participant_counts.tsv"
  echo "  output/tables/phenotype_missingness_summary.tsv"
  echo "  output/tables/pheno_with_pcs.txt"
  echo "  output/tables/sum_stats.txt"
  echo "  output/tables/manhattan_input.txt"
  echo "  output/tables/sum_stats_AF_lt_0.01.txt"
  echo "  output/tables/sum_stats_AF_ge_0.01.txt"
  echo "  output/tables/sum_stats_AF_lt_0.01_CHR_BP_SNP_PVAL.txt"
  echo "  output/tables/sum_stats_AF_ge_0.01_CHR_BP_SNP_PVAL.txt"
  echo "  output/plots/Rect_Manhtn.common_plot_input.jpg"
  echo "  output/plots/common_plot_input_qq_lambda.jpg"
  echo "  output/plots/common_plot_input_circular.jpg"
  echo "  output/plots/common_plot_input_density.jpg"
  echo "  output/plots/Rect_Manhtn.rare_plot_input.jpg"
  echo "  output/plots/rare_plot_input_qq_lambda.jpg"
  echo "  output/plots/rare_plot_input_circular.jpg"
  echo "  output/plots/rare_plot_input_density.jpg"
  echo "  output/plots/common_plots_skipped.txt"
  echo "  output/plots/rare_plots_skipped.txt"
  echo "  config/run_manifest.txt"
  echo "  logs/pipeline.log"
  echo "  qc_summary.tsv"
  echo "  diagnostics/covariate_sanity_report.tsv"
  echo "  diagnostics/null_model_diagnostics.tsv"
  echo
  echo "Model configuration:"
  echo "  matched samples         : ${N_COMMON}"
  echo "  use sparse GRM          : ${USE_SPARSE_GRM}"
  echo "  sparse GRM threshold    : ${SPARSE_GRM_MIN_SAMPLES}"
  echo "  Firth p cutoff          : ${P_CUTOFF_FOR_FIRTH}"
  echo "  association MAF filter  : ${ASSOCIATION_MAF_THRESHOLD}"
  echo "  rare/common split AF    : ${COMMON_RARE_SPLIT_THRESHOLD}"
  echo "  cleanup on success      : ${CLEANUP_ON_SUCCESS}"
  echo
  echo "Sample counts:"
  echo "  matched.assoc fam rows   : ${N_ASSOC}"
  echo "  matched.kinship fam rows : ${N_KINSHIP}"
  echo
  echo "Plot input counts:"
  echo "  common variants          : $(count_data_rows "${MANHATTAN_INPUT_COMMON}")"
  echo "  rare variants            : $(count_data_rows "${MANHATTAN_INPUT_RARE}")"
} > output/tables/run_summary.txt

echo "Step 15. Cleanup temporary/intermediate files after successful completion..."
CLEANUP_REPORT="output/tables/cleanup_report.txt"

if [ "${CLEANUP_ON_SUCCESS}" = "true" ]; then
  {
    echo "Cleanup performed on: $(date '+%Y-%m-%d %H:%M:%S')"
    echo
    echo "Removed temporary and intermediate files:"
  } > "${CLEANUP_REPORT}"

  cleanup_path() {
    local path=$1
    if [ -e "${path}" ]; then
      rm -rf "${path}"
      echo "  ${path}" >> "${CLEANUP_REPORT}"
    fi
  }

  cleanup_glob() {
    local pattern=$1
    for item in ${pattern}; do
      if [ -e "${item}" ]; then
        rm -rf "${item}"
        echo "  ${item}" >> "${CLEANUP_REPORT}"
      fi
    done
    return 0
  }

  cleanup_glob "QCed.assoc.*"
  cleanup_glob "QCed.kinship.*"
  cleanup_glob "matched.assoc.*"
  cleanup_glob "matched.kinship.*"

  cleanup_path "intermediate/pheno.cleaned.all.txt"
  cleanup_path "intermediate/pheno.analysis_ready.txt"
  cleanup_glob "intermediate/mypc.*"

  cleanup_path "output/plots/common_plot_input.txt"
  cleanup_path "output/plots/rare_plot_input.txt"

  if [ -d "sparseGRM" ]; then
    cleanup_glob "sparseGRM/*"
  fi

  echo >> "${CLEANUP_REPORT}"
  echo "Kept final deliverables in output/tables, output/plots, logs, config, and diagnostics." >> "${CLEANUP_REPORT}"
else
  {
    echo "Cleanup skipped because CLEANUP_ON_SUCCESS=${CLEANUP_ON_SUCCESS}"
  } > "${CLEANUP_REPORT}"
fi

echo "Pipeline finished successfully."
echo "See output/tables/run_summary.txt for a compact summary."
echo "See output/tables/cleanup_report.txt for cleanup details."