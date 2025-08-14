#!/bin/bash

WORKING_DIR=$1
KINSHIP_INPUT=$2
DOSAGE_INPUT=$3
PHENO_INPUT=$4
MODEL=$5
TYPE=$6


# Add pipeline bin folder to PATH
export PATH="$(pwd)/bin:$PATH"

# Run your pipeline
bash submit_all.sh "$WORKING_DIR" "$KINSHIP_INPUT" "$DOSAGE_INPUT" "$PHENO_INPUT" "$MODEL" "$TYPE"

# organize
#mkdir output
cd output
mkdir plots tables
cd ../
mv *jpg ./output/plots
mv manhattan_plot_input_ready.txt pheno_with_pcs.txt  sum_stats.txt ./output/tables
rm -r logs tmp
rm my*
