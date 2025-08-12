#!/bin/bash


INPUT_GENO=$1
INPUT_PHENO=$2
COVAR_LIST=$3
QCOVAR_LIST=$4
DIAG_INPUT=$5
TYPE=$6


CHR=$5  # Single chromosome to process (e.g. 1 to 22)

export PATH="$(pwd)/bin:$PATH"

# Start timer
START_TIME=$(date +%s)

# Call the analysis script
bash saige.sh "$INPUT_GENO" "$INPUT_PHENO" "$COVAR_LIST" "$QCOVAR_LIST" "$DIAG_INPUT" "$TYPE"

# Stop timer
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

# Print and log runtime
echo "Runtime: $RUNTIME seconds"
echo "$(date): ${GENO_INPUT}, ${PHENO_INPUT} -> ${RUNTIME} seconds" >> runtime_log.txt

# organize
mkdir output
cd output
mkdir plots tables
cd ../
mv *jpg ./output/plots
mv ./saige_output/manhattan_input.txt pheno_with_pcs.txt  sum_stats.txt ./output/tables
#mv ./saige_output/manhattan_input.txt ./output/tables/
rm -r logs saige_output/ sparseGRM/ summary_plots/
rm mypc.*



