#!/bin/bash

WORKING_DIR=$1
INPUT_PATH=$2
MODEL=$4


# Add pipeline bin folder to PATH
export PATH="$(pwd)/bin:$PATH"

# Run your pipeline
bash submit_all.sh "$WORKING_DIR" "$INPUT_PATH" "$INPUT_PATH" "$MODEL"

# example how to run it: 
#   ./run_pipeline_gmmat.sh ./ /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/PLINK_FILES/PRADI/ /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/PLINK_FILES/PRADI/ "ADRD ~ AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5"