#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="${SAGA_GMMAT_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BIN_DIR="$SCRIPT_DIR/../bin"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/tmp"

RAW_PHENO="$PHENO_INPUT"
CLEAN_PHENO="$OUTPUT_DIR/pheno_cleaned.txt"

ASSOC_BFILE="$OUTPUT_DIR/QCed.assoc"
KINSHIP_BFILE="$OUTPUT_DIR/QCed.kinship"

ASSOC_COMMON_PREFIX="$OUTPUT_DIR/QCed.assoc.common"
KINSHIP_COMMON_PREFIX="$OUTPUT_DIR/QCed.kinship.common"

COMMON_IDS="$OUTPUT_DIR/common_ids_in_order.txt"
PHENO_FILTERED="$OUTPUT_DIR/pheno_filtered_common.txt"
PHENO_GEMMA="$OUTPUT_DIR/pheno_fid_kinship.txt"
PHENO_WITH_PCS="$OUTPUT_DIR/pheno_with_pcs.txt"

ASSOC_IDS="$OUTPUT_DIR/qced_assoc_ids.tmp"
KINSHIP_IDS="$OUTPUT_DIR/qced_kinship_ids.tmp"
PHENO_IDS="$OUTPUT_DIR/pheno_ids.tmp"

KINSHIP_OUT_PREFIX="mykinship"
PCA_OUT="$OUTPUT_DIR/mypc"

echo ">> Running QC on genotype data"
bash "$SCRIPT_DIR/qc.sh" "$GENO_INPUT" "$GENO_INPUT" "$OUTPUT_DIR"

if [ ! -f "$ASSOC_BFILE.fam" ]; then
    echo "ERROR: expected file not found: $ASSOC_BFILE.fam"
    exit 1
fi

if [ ! -f "$KINSHIP_BFILE.fam" ]; then
    echo "ERROR: expected file not found: $KINSHIP_BFILE.fam"
    exit 1
fi

echo ">> Cleaning phenotype file and replacing missing values with NA"
python3 - "$RAW_PHENO" "$CLEAN_PHENO" <<'PY'
import csv
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r", newline="") as fin, open(outfile, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

    for row_idx, row in enumerate(reader):
        cleaned = []
        for value in row:
            if value is None:
                cleaned.append("NA")
            else:
                stripped = value.strip()
                if stripped == "" or stripped == ".":
                    cleaned.append("NA")
                else:
                    cleaned.append(stripped)
        writer.writerow(cleaned)
PY

echo ">> Extracting sample IDs from QCed.assoc.fam"
awk 'BEGIN{OFS="\t"} {print $1, $2}' "$ASSOC_BFILE.fam" > "$ASSOC_IDS"

echo ">> Extracting sample IDs from QCed.kinship.fam"
awk 'BEGIN{OFS="\t"} {print $1, $2}' "$KINSHIP_BFILE.fam" > "$KINSHIP_IDS"

echo ">> Extracting sample IDs from cleaned phenotype file"
awk '
BEGIN {
    FS = OFS = "\t"
}
NR == 1 {
    fid_col = 0
    iid_col = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "FID") {
            fid_col = i
        } else if ($i == "IID") {
            iid_col = i
        }
    }
    if (fid_col == 0 || iid_col == 0) {
        print "ERROR: phenotype file must contain FID and IID columns" > "/dev/stderr"
        exit 1
    }
    next
}
{
    print $fid_col, $iid_col
}
' "$CLEAN_PHENO" > "$PHENO_IDS"

echo ">> Building common sample intersection across QCed.assoc, QCed.kinship, and cleaned phenotype"
awk '
BEGIN {
    FS = OFS = "\t"
}
FNR == NR {
    assoc[$1 OFS $2] = 1
    next
}
ARGIND == 2 {
    kinship[$1 OFS $2] = 1
    next
}
ARGIND == 3 {
    pheno[$1 OFS $2] = 1
    next
}
ARGIND == 4 {
    key = $1 OFS $2
    if ((key in assoc) && (key in kinship) && (key in pheno)) {
        print $1, $2
    }
}
' "$ASSOC_IDS" "$KINSHIP_IDS" "$PHENO_IDS" "$KINSHIP_IDS" > "$COMMON_IDS"

N_ASSOC=$(wc -l < "$ASSOC_IDS")
N_KINSHIP=$(wc -l < "$KINSHIP_IDS")
N_PHENO_IDS=$(wc -l < "$PHENO_IDS")
N_COMMON=$(wc -l < "$COMMON_IDS")

echo ">> Samples in QCed.assoc.fam: $N_ASSOC"
echo ">> Samples in QCed.kinship.fam: $N_KINSHIP"
echo ">> Samples in cleaned phenotype file: $N_PHENO_IDS"
echo ">> Samples in common intersection: $N_COMMON"

if [ "$N_COMMON" -le 0 ]; then
    echo "ERROR: no common samples found across QCed.assoc, QCed.kinship, and phenotype"
    exit 1
fi

echo ">> Subsetting QCed.assoc to common sample set"
"$BIN_DIR/plink" \
    --bfile "$ASSOC_BFILE" \
    --keep "$COMMON_IDS" \
    --make-bed \
    --out "$ASSOC_COMMON_PREFIX"

echo ">> Subsetting QCed.kinship to common sample set"
"$BIN_DIR/plink" \
    --bfile "$KINSHIP_BFILE" \
    --keep "$COMMON_IDS" \
    --make-bed \
    --out "$KINSHIP_COMMON_PREFIX"

for ext in bed bim fam; do
    mv -f "$ASSOC_COMMON_PREFIX.$ext" "$ASSOC_BFILE.$ext"
    mv -f "$KINSHIP_COMMON_PREFIX.$ext" "$KINSHIP_BFILE.$ext"
done

rm -f "$ASSOC_COMMON_PREFIX.log" "$ASSOC_COMMON_PREFIX.nosex"
rm -f "$KINSHIP_COMMON_PREFIX.log" "$KINSHIP_COMMON_PREFIX.nosex"

echo ">> Filtering cleaned phenotype file to common sample set and preserving genotype order"
awk '
BEGIN {
    FS = OFS = "\t"
}
NR == FNR {
    keep_order[NR] = $1 OFS $2
    keep_seen[$1 OFS $2] = 1
    n_keep = NR
    next
}
FNR == 1 {
    fid_col = 0
    iid_col = 0
    pheno_col = 0

    for (i = 1; i <= NF; i++) {
        if ($i == "FID") {
            fid_col = i
        } else if ($i == "IID") {
            iid_col = i
        } else if ($i == "PHENO") {
            pheno_col = i
        }
    }

    if (fid_col == 0 || iid_col == 0 || pheno_col == 0) {
        print "ERROR: phenotype file must contain FID, IID, and PHENO columns" > "/dev/stderr"
        exit 1
    }

    header_line = $1
    for (i = 2; i <= NF; i++) {
        header_line = header_line OFS $i
    }
    next
}
{
    key = $fid_col OFS $iid_col
    if (key in keep_seen) {
        row_by_id[key] = $0
    }
}
END {
    print header_line

    missing_count = 0
    for (i = 1; i <= n_keep; i++) {
        key = keep_order[i]
        if (key in row_by_id) {
            print row_by_id[key]
        } else {
            print "ERROR: sample in common ID list not found in phenotype file: " key > "/dev/stderr"
            missing_count++
        }
    }

    if (missing_count > 0) {
        exit 1
    }
}
' "$COMMON_IDS" "$CLEAN_PHENO" > "$PHENO_FILTERED"

echo ">> Validating filtered phenotype column consistency"
awk -F'\t' '
NR == 1 {
    expected = NF
    next
}
NF != expected {
    print "ERROR: filtered phenotype file has inconsistent number of columns at line " NR ": expected " expected ", found " NF > "/dev/stderr"
    exit 1
}
' "$PHENO_FILTERED"

N_PHENO_FILTERED=$(awk 'END{print NR-1}' "$PHENO_FILTERED")
echo ">> Samples in filtered phenotype file: $N_PHENO_FILTERED"

if [ "$N_PHENO_FILTERED" -ne "$N_COMMON" ]; then
    echo "ERROR: filtered phenotype sample count does not match common intersection"
    exit 1
fi

echo ">> Creating GEMMA phenotype file"
awk '
BEGIN {
    FS = OFS = "\t"
}
NR == 1 {
    pheno_col = 0
    for (i = 1; i <= NF; i++) {
        if ($i == "PHENO") {
            pheno_col = i
        }
    }
    if (pheno_col == 0) {
        print "ERROR: filtered phenotype file must contain PHENO column" > "/dev/stderr"
        exit 1
    }
    next
}
{
    print $pheno_col
}
' "$PHENO_FILTERED" > "$PHENO_GEMMA"

echo ">> Computing kinship matrix with GEMMA"
(
    cd "$OUTPUT_DIR"
    "$BIN_DIR/gemma-0.98.5" \
        -bfile "QCed.kinship" \
        -gk 1 \
        -p "$(basename "$PHENO_GEMMA")" \
        -maf 0.05 \
        -o "$KINSHIP_OUT_PREFIX"
)

if [ -d "$OUTPUT_DIR/output" ]; then
    shopt -s nullglob
    gemma_files=("$OUTPUT_DIR"/output/*)
    if [ "${#gemma_files[@]}" -gt 0 ]; then
        mv "${gemma_files[@]}" "$OUTPUT_DIR/"
    fi
    shopt -u nullglob
    rmdir "$OUTPUT_DIR/output" 2>/dev/null || true
fi

echo ">> Computing top 10 PCs with PLINK"
"$BIN_DIR/plink" --bfile "$KINSHIP_BFILE" --pca 10 --out "$PCA_OUT"

if [ ! -f "${PCA_OUT}.eigenvec" ]; then
    echo "ERROR: expected PCA output not found: ${PCA_OUT}.eigenvec"
    exit 1
fi

echo ">> Merging PCs with filtered phenotype"
awk '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
}
NR == FNR {
    if (NR == 1) {
        pheno_header = $0
        next
    }
    key = $1 OFS $2
    pheno_row[key] = $0
    next
}
FNR == 1 {
    print pheno_header, "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
}
{
    key = $1 OFS $2
    if (!(key in pheno_row)) {
        print "ERROR: PCA sample not found in filtered phenotype file: " key > "/dev/stderr"
        exit 1
    }
    print pheno_row[key], $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
}
END {
    if (pheno_header == "") {
        print "ERROR: filtered phenotype file header was not captured" > "/dev/stderr"
        exit 1
    }
}
' "$PHENO_FILTERED" "${PCA_OUT}.eigenvec" > "$PHENO_WITH_PCS"

echo ">> Validating pheno_with_pcs column consistency"
awk -F'\t' '
NR == 1 {
    expected = NF
    next
}
NF != expected {
    print "ERROR: pheno_with_pcs.txt has inconsistent number of columns at line " NR ": expected " expected ", found " NF > "/dev/stderr"
    exit 1
}
' "$PHENO_WITH_PCS"

N_PCS=$(awk 'END{print NR-1}' "$PHENO_WITH_PCS")
echo ">> Samples in phenotype with PCs: $N_PCS"

if [ "$N_PCS" -ne "$N_COMMON" ]; then
    echo "ERROR: pheno_with_pcs sample count does not match common intersection"
    exit 1
fi

echo ">> Step 1 completed successfully"