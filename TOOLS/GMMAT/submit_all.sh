#!/bin/bash

# Check arguments
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <working_dir> <BFILE_PATH_MAF0_05_R2_099> <DOSAGE_DIR_MAF0_01_R2_080> <MODEL_TYPE>"
  exit 1
fi

# Arguments
WORKDIR="$1"
BFILE_PATH_MAF0_05_R2_099="$2"
DOSAGE_DIR_MAF0_01_R2_080="$3"
MODEL_TYPE="$4"

# Ensure logs directory exists
mkdir -p "$WORKDIR/logs"

# Submit step 1
job1=$(sbatch --parsable --export=ALL,WORKDIR="$WORKDIR",BFILE_PATH_MAF0_05_R2_099="$BFILE_PATH_MAF0_05_R2_099" step1_run_kinship_pca.sh)
echo "Submitted step1_run_kinship_pca.sh as job $job1"

# Submit step 2
#MODEL_TYPE="model1"  # or model2
job2=$(sbatch --parsable \
  --dependency=afterok:$job1 \
  --export=ALL,WORKDIR="$WORKDIR",BFILE_PATH_MAF0_05_R2_099="$BFILE_PATH_MAF0_05_R2_099",DOSAGE_DIR_MAF0_01_R2_080="$DOSAGE_DIR_MAF0_01_R2_080",MODEL_TYPE="$MODEL_TYPE" \
  step2_run_full_gwas_pipeline.slurm)

# Submit step 3
job3=$(sbatch --parsable --dependency=afterok:$job2 --export=ALL,WORKDIR="$WORKDIR" step3_postprocess_gwas.sh)
echo "Submitted step3_postprocess_gwas.sh as job $job3 (after $job2)"


# How to run it:
#   bash submit_all.sh ./ /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/PLINK_FILES/ADC_MEX_CH/ /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/PLINK_FILES/ADC_MEX_CH/ model1