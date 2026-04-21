#!/usr/bin/env bash
#===========================================#
#         GWAS PIPELINE (PLINK + R)         #
#===========================================#

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ $# -lt 6 ]; then
  echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <COVAR_LIST> <DIAGNOST_VAR> <TRAIT_TYPE> <OUTPUT_DIR>"
  exit 1
fi

GENO_INPUT="$1"
INPUT_PHENO="$2"
COVAR_LIST="$3"
DIAGNOST_VAR="$4"
TRAIT_TYPE="$5"
OUTPUT_DIR="$6"

resolve_plink() {
  if [ -n "${PLINK_BIN:-}" ] && [ -x "${PLINK_BIN}" ]; then
    echo "${PLINK_BIN}"
    return 0
  fi

  if command -v plink >/dev/null 2>&1; then
    command -v plink
    return 0
  fi

  if [ -x "${SCRIPT_DIR}/../bin/plink" ]; then
    echo "${SCRIPT_DIR}/../bin/plink"
    return 0
  fi

  echo "ERROR: Could not find PLINK. Set PLINK_BIN or add plink to PATH." >&2
  exit 1
}

PLINK="$(resolve_plink)"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "ERROR: Rscript not found in PATH."
  exit 1
fi

plink_prefix_exists() {
  local prefix="$1"
  [ -f "${prefix}.bed" ] && [ -f "${prefix}.bim" ] && [ -f "${prefix}.fam" ]
}

require_file() {
  local f="$1"
  if [ ! -f "$f" ]; then
    echo "ERROR: Required file not found: $f" >&2
    exit 1
  fi
}

require_nonempty_data_file() {
  local f="$1"
  local label="$2"

  require_file "$f"

  local nlines
  nlines=$(wc -l < "$f")

  if [ "$nlines" -le 1 ]; then
    echo "ERROR: ${label} has no data rows: $f" >&2
    echo "First lines of ${f}:" >&2
    head -5 "$f" >&2 || true
    exit 1
  fi
}

print_count() {
  local label="$1"
  local f="$2"
  if [ -f "$f" ]; then
    echo "${label}: $(wc -l < "$f") lines -> $f"
  else
    echo "${label}: missing -> $f"
  fi
}

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "All outputs will be stored in: $OUTPUT_DIR"
echo "Using PLINK: ${PLINK}"
echo "Using Rscript: $(command -v Rscript)"

if [[ "$TRAIT_TYPE" == "binary" ]]; then
  TEST_ARGS=(--logistic hide-covar --allow-no-sex --1 --freq)
  GWAS_FILE="assoc.assoc.logistic"
  EFFECT_HEADER="OR"
elif [[ "$TRAIT_TYPE" == "quantitative" ]]; then
  TEST_ARGS=(--linear hide-covar --allow-no-sex --freq)
  GWAS_FILE="assoc.assoc.linear"
  EFFECT_HEADER="BETA"
else
  echo "Invalid TRAIT_TYPE: must be 'binary' or 'quantitative'"
  exit 1
fi

echo "Running association analysis for ${TRAIT_TYPE} trait"
echo "Input genotype dataset:"
echo "${GENO_INPUT}"

#===========================================#
#        0. QUALITY CONTROL (QC)            #
#===========================================#
if plink_prefix_exists "QCed.assoc" && plink_prefix_exists "QCed.kinship"; then
  echo "Skipping QC (QCed.assoc and QCed.kinship already exist)"
else
  "${SCRIPT_DIR}/qc.sh" "${GENO_INPUT}" "${GENO_INPUT}"
fi

require_file "QCed.assoc.bed"
require_file "QCed.assoc.bim"
require_file "QCed.assoc.fam"
require_file "QCed.kinship.bed"
require_file "QCed.kinship.bim"
require_file "QCed.kinship.fam"

echo "QC complete"
print_count "QCed.assoc.fam" "QCed.assoc.fam"
print_count "QCed.kinship.fam" "QCed.kinship.fam"

GENO_INPUT="./QCed.assoc"
PCA_INPUT="./QCed.kinship"

#===========================================#
#                1. PCA                     #
#===========================================#
PCA_OUT="./mypc"

if [ -f "${PCA_OUT}.eigenvec" ]; then
  echo "Skipping PCA (${PCA_OUT}.eigenvec already exists)"
else
  echo "PCA input:"
  echo "${PCA_INPUT}"
  echo "Current directory: $(pwd)"
  "${PLINK}" --bfile "${PCA_INPUT}" --pca 10 --out "${PCA_OUT}"
fi

require_file "${PCA_OUT}.eigenvec"
require_file "${PCA_OUT}.eigenval"

echo "PCA complete"
print_count "PCA eigenvec" "${PCA_OUT}.eigenvec"

#===========================================#
#    2. HARMONIZE PHENO + PCS WITH R        #
#===========================================#
PHENO_MATCHED="./pheno_matched.txt"
PHENO_WITH_PCS="./pheno_with_pcs.txt"

if [ -f "${PHENO_MATCHED}" ] && [ -f "${PHENO_WITH_PCS}" ]; then
  echo "Skipping phenotype/PCA merge (${PHENO_MATCHED} and ${PHENO_WITH_PCS} already exist)"
else
  Rscript "${SCRIPT_DIR}/merge_pheno_pcs.R" \
    "${INPUT_PHENO}" \
    "${GENO_INPUT}.fam" \
    "${PCA_INPUT}.fam" \
    "${PCA_OUT}.eigenvec" \
    "${PHENO_MATCHED}" \
    "${PHENO_WITH_PCS}"
fi

require_nonempty_data_file "${PHENO_MATCHED}" "Matched phenotype file"
require_nonempty_data_file "${PHENO_WITH_PCS}" "Phenotype with PCs file"

echo "Phenotype/PCA merge complete"
print_count "pheno_matched.txt" "${PHENO_MATCHED}"
print_count "pheno_with_pcs.txt" "${PHENO_WITH_PCS}"

#===========================================#
#         3. ASSOCIATION TESTING            #
#===========================================#
if [ -f "./${GWAS_FILE}" ] && [ -f "./assoc.frq" ]; then
  echo "Skipping association step (GWAS output already exists)"
else
  "${PLINK}" \
    --bfile "${GENO_INPUT}" \
    --pheno "${PHENO_WITH_PCS}" \
    --pheno-name "${DIAGNOST_VAR}" \
    --covar "${PHENO_WITH_PCS}" \
    --covar-name "${COVAR_LIST}" \
    "${TEST_ARGS[@]}" \
    --out "./assoc"
fi

require_nonempty_data_file "./${GWAS_FILE}" "Association result file"
require_nonempty_data_file "./assoc.frq" "Allele frequency file"

echo "Association complete"
print_count "${GWAS_FILE}" "./${GWAS_FILE}"
print_count "assoc.frq" "./assoc.frq"

#===========================================#
#      4. MERGE FREQUENCIES + GWAS          #
#===========================================#
FREQ_FILE="./assoc.frq"
GWAS_RESULT_FILE="./${GWAS_FILE}"
MERGED_FILE="./assoc.merged.txt"

rm -f "${MERGED_FILE}" "./sum_stats.txt" "./manhattan_input.txt"

echo "Running merge step..."
ls -lh "${FREQ_FILE}" "${GWAS_RESULT_FILE}"

awk -v effect_header="${EFFECT_HEADER}" '
  BEGIN { FS = "[ \t]+"; OFS = "\t" }

  NR == FNR {
    sub(/^[ \t]+/, "", $0)
    if (FNR == 1) next
    if (NF < 6) next
    maf[$2] = $5
    nchrobs[$2] = $6
    next
  }

  {
    sub(/^[ \t]+/, "", $0)

    if (FNR == 1) {
      print "CHR", "SNP", "BP", "A1", "TEST", "NMISS", effect_header, "STAT", "PVAL", "MAF", "NCHROBS"
      next
    }

    if (NF < 9) next

    if ($2 in maf && $2 in nchrobs) {
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, maf[$2], nchrobs[$2]
    }
  }
' "${FREQ_FILE}" "${GWAS_RESULT_FILE}" > "${MERGED_FILE}"

require_nonempty_data_file "${MERGED_FILE}" "Merged GWAS+frequency file"

echo "Merge completed"
print_count "assoc.merged.txt" "${MERGED_FILE}"
head -5 "${MERGED_FILE}"

#===========================================#
#        5. POST-PROCESSING RESULTS         #
#        FILTER OUT INVALID / NA PVAL       #
#===========================================#
awk '
  BEGIN { FS = "[ \t]+"; OFS = "\t" }

  FNR == 1 {
    print $0
    next
  }

  $9 != "NA" && $9 != "" && $9 != "nan" && \
  $10 != "NA" && $10 != "" && $10 != "nan" && \
  ($9 + 0) > 0 && ($9 + 0) <= 1 && \
  ($10 + 0) >= 0.01 && ($10 + 0) < 0.99 {
    print
  }
' "${MERGED_FILE}" > "./sum_stats.txt"

require_nonempty_data_file "./sum_stats.txt" "Filtered summary statistics file"

echo "Post-processing completed"
print_count "sum_stats.txt" "./sum_stats.txt"
head -5 "./sum_stats.txt"

#===========================================#
#     6. MANHATTAN INPUT + SANITY CHECKS    #
#===========================================#
awk '
  BEGIN { FS = "[ \t]+"; OFS = "\t" }

  NR == 1 {
    print "CHR", "BP", "SNP", "PVAL"
    next
  }

  $9 != "NA" && $9 != "" && ($9 + 0) > 0 && ($9 + 0) <= 1 {
    print $1, $3, $2, $9
  }
' "./sum_stats.txt" > "./manhattan_input.txt"

require_nonempty_data_file "./manhattan_input.txt" "Manhattan input file"

echo "Manhattan input created"
print_count "manhattan_input.txt" "./manhattan_input.txt"
head -5 "./manhattan_input.txt"

awk '
  BEGIN { FS = "[ \t]+" }
  NR > 1 {
    if (NF != 4) {
      print "ERROR: manhattan_input.txt has a line with " NF " columns instead of 4 at line " NR > "/dev/stderr"
      exit 1
    }
    if ($4 == "" || $4 == "NA" || ($4 + 0) <= 0 || ($4 + 0) > 1) {
      print "ERROR: manhattan_input.txt has invalid PVAL at line " NR ": " $4 > "/dev/stderr"
      exit 1
    }
  }
' "./manhattan_input.txt"

#===========================================#
#           7. MANHATTAN + QQ PLOTS         #
#===========================================#
Rscript "${SCRIPT_DIR}/create_manhattan.R" manhattan_input.txt
Rscript "${SCRIPT_DIR}/create_qq.plot.R" manhattan_input.txt
Rscript "${SCRIPT_DIR}/create_circular_manhattan.R" manhattan_input.txt
Rscript "${SCRIPT_DIR}/create_density_plot.R" manhattan_input.txt

echo "Plotting completed"

#===========================================#
#               8. ORGANIZING               #
#===========================================#
mkdir -p output/plots output/tables

if [ -f "Cir_Manhtn.manhattan_input_circular.jpg" ]; then
  mv -f "Cir_Manhtn.manhattan_input_circular.jpg" "output/plots/circular_manhattan_plot.jpg"
fi

if [ -f "Rect_Manhtn.manhattan_input_manhattan_highlight.jpg" ]; then
  mv -f "Rect_Manhtn.manhattan_input_manhattan_highlight.jpg" "output/plots/rectangular_manhattan_plot.jpg"
fi

if [ -f "manhattan_input_qq_lambda.jpg" ]; then
  mv -f "manhattan_input_qq_lambda.jpg" "output/plots/Q_Q_plot.jpg"
fi

if [ -f "Marker_Density.manhattan_input_density.jpg" ]; then
  mv -f "Marker_Density.manhattan_input_density.jpg" "output/plots/SNP_density_plot.jpg"
fi

if [ -f "manhattan_input.txt" ]; then
  cp -f "manhattan_input.txt" "output/tables/manhattan_input.txt"
fi

if [ -f "pheno_with_pcs.txt" ]; then
  cp -f "pheno_with_pcs.txt" "output/tables/pheno_with_pcs.txt"
fi

if [ -f "pheno_matched.txt" ]; then
  cp -f "pheno_matched.txt" "output/tables/pheno_matched.txt"
fi

if [ -f "sum_stats.txt" ]; then
  cp -f "sum_stats.txt" "output/tables/sum_stats.txt"
fi

echo "Processing complete. Results are in: ${OUTPUT_DIR}"
pwd
mv Rect_Manhtn.manhattan_input.jpg output/plots/
