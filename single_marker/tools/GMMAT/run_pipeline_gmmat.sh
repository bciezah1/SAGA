#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <GENO_INPUT> <PHENO_INPUT> <MODEL> <TYPE> <OUTPUT_DIR>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

GENO_INPUT="$1"
PHENO_INPUT="$2"
MODEL="$3"
TYPE="$4"
OUTPUT_DIR="$5"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/tmp"

export GENO_INPUT
export PHENO_INPUT
export MODEL_FORMULA="$MODEL"
export TYPE
export OUTPUT_DIR
export SAGA_GMMAT_SCRIPT_DIR="$SCRIPT_DIR"

bash "$SCRIPT_DIR/submit_all.sh" "$GENO_INPUT" "$PHENO_INPUT" "$MODEL" "$TYPE" "$OUTPUT_DIR"