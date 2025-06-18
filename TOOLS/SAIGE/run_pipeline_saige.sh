#!/bin/bash

# Check if path argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <input_path>"
  exit 1
fi

INPUT_PATH=$1
MODEL=$2

# Load conda
#source ~/miniconda3/etc/profile.d/conda.sh

# Create env if it doesn't exist
#if ! conda env list | grep -q "^pipeline_env"; then
#    conda env create -f environment.yml
#fi

#conda activate pipeline_env

# Add pipeline bin folder to PATH
export PATH="$(pwd)/bin:$PATH"

# Run your pipeline
sbatch saige.slurm "$INPUT_PATH" "$MODEL"

# example - how to run it:
#./run_pipeline_saige.sh ../../toy_data/ model1

#conda deactivate
