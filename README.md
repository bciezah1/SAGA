# EeasyGWAS 

        ## GMMAT

        run command:
                - bash submit_all.sh ./ ./../toy_data/ ../toy_data/ model1

        ## SAIGE

                - sbatch saige.slurm ../toy_data/ model1

This repository contains a SLURM-based pipeline for conducting single-marker genome-wide association studies (GWAS) using [GMMAT](https://github.com/hanchenphd/GMMAT), [PLINK](https://www.cog-genomics.org/plink/), and [GEMMA](https://github.com/genetics-statistics/GEMMA).



## 🧾 Required Inputs

| Parameter | Description |
|----------|-------------|
| `WORKDIR` | Full path to working directory |
| `BFILE_PATH_MAF0_05_R2_099` | Directory containing merged PLINK file: `chr1_chr22.r2.99.maf0.05.{bed,bim,fam}` |
| `DOSAGE_DIR_MAF0_01_R2_080` | Directory containing per-chromosome dosage files: `chr1.r2.80.maf0.01`, ..., `chr22.r2.80.maf0.01` |
| `MODEL_TYPE` | Either `model1` or `model2` (see [Model Definitions](#model-definitions)) |

---

## 🚀 How to Run

From your terminal:

```bash

Generic command line:
bash submit_all.sh <WORKDIR> <BFILE_PATH_MAF0_05_R2_099> <DOSAGE_DIR_MAF0_01_R2_080> <MODEL_TYPE>

Quick start:
bash submit_all.sh ./ ./toy_data/ ./toy_data/ model1

```

## 🧠 Model Definitions
model1
ADRD ~ AGE + as.factor(SEX) + PC1 + PC2 + PC3 + PC4 + PC5

model2
ADRD ~ AGE + as.factor(SEX) + PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(APOEe4)

## ⚙️ Pipeline Steps
#### Step 1: Kinship & PCA (step1_run_kinship_pca.sh)

Steps: 

    - Extracts FID/IID/ADRD from phenotype file
    - Computes kinship matrix using GEMMA
    - Computes PCs using PLINK
    - Merges PCs with phenotype data

Outputs:

    - pheno_with_pcs.txt
    - mykinship.cXX.txt
    - mypc.eigenvec, mypc.eigenval

#### Step 2: GWAS (step2_run_full_gwas_pipeline.slurm)

Steps:

    - Runs a SLURM job array over chromosomes 1–22:
    - Builds GLMM model with kinship matrix
    - Performs score test using glmm.score()

Outputs (per chromosome):

    - raw summary statistics

#### Step 3: Postprocessing (step3_postprocess_gwas.sh)

Steps:

    - Concatenates all output
    - Calculates BETA and SE
    - Filters by allele frequency (0.01 ≤ AF < 0.99)
    - Generates Manhattan and QQ plots

Outputs:

    - clean statistics
    - manhattan plot
    - qqplot
