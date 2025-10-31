#!/bin/bash

#WORKING_DIR=$1
GENO_INPUT=$1
PHENO_INPUT=$2
MODEL=$3
TYPE=$4
OUTPUT_DIR=$5

# Create output folder
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/tmp"

# Export variables for sub-scripts
#export WORKDIR="$WORKING_DIR"
export GENO_INPUT
export PHENO_INPUT
export MODEL_FORMULA="$MODEL"
export TYPE
export OUTPUT_DIR

# Add bin folder to PATH
export PATH="$(pwd)/bin:$PATH"

# Run main pipeline
bash submit_all.sh "$GENO_INPUT" "$PHENO_INPUT" "$MODEL" "$TYPE" "$OUTPUT_DIR"
