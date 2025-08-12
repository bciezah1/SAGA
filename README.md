# SAGA

## 🧬 Summary

**SAGA** is a collection of streamlined pipelines for performing single-marker GWAS using **PLINK**, **GMMAT** and **SAIGE**. It is designed for scientists with limited programming experience who want to analyze genotyped or imputed genetic data with minimal setup.

---

## 📥 Required Inputs

SAGE requires **three main input files**. To minimize errors, we recommend renaming your files to match the naming conventions shown in the `toy_data/` folder.

### Input Files:

- **Kinship input (PLINK format):** High-quality imputed genotype data for kinship estimation.  
  _Recommended: `R² ≥ 0.99`, `MAF ≥ 0.05`_

  expected label file:
  
      input_kinship.bed
      input_kinship.bim
      input_kinship.fam

- **Association input (PLINK format):** Dosage or genotype data for association testing.  
  _Recommended: `R² ≥ 0.80`, `MAF ≥ 0.01`_

  expected label files:

      input_dosage.bed
      input_dosage.bim
      input_dosage.fam

- **Phenotype file:** Must include the following columns:  
  `FID`, `IID`, `PHENO`, `SEX`, `AGE`, `APOEe4`

  expected pheno file format:

      FID     IID     SEX     AGE     PHENO
      FAM001  IND001  0       84      0
      FAM002  IND002  0       85      1
      FAM003  IND003  1       72      0

      ...


> ⚠️ **Warning:** The number of participants and their order **must match exactly** between the PLINK files and the phenotype file.

> ⚠️ **Important:** Refer to the `toy_data/` folder to verify correct formatting and file naming.

> ⚠️ **Important:** You will need to have installed on your linux system:
        R (Version 4.2.2) which should include the next libraries: GMMAT, CMplot, qqman

>  ⚠️ **Important:** Open permit to executables
       Open permits for all files on the bin directory
       Open permits for all files inside each pipeline folder GMMAT, SAIGE, and PLINK
       Tip: To open permits, just go inside the folder mentioned, and run the next command: chmod +x *
         
---

## 🚀 Pipelines

SAGA includes three pipelines: The first one using **PLINK**, the second one using **GMMAT**, and third one using **SAIGE**. We recommend running the provided toy data first to verify that everything is working properly before applying the pipeline to your own dataset.

### 🔧 PLINK Pipeline

**Run with:**

```bash

1. Get inside the PLINK folder.
2. Run the command:
      ./run_pipeline_plink.sh ../../../toy_data/genotype ../../../toy_data/pheno_continue.txt  SEX,AGE PHENO quantitative


# Explanation

./run_pipeline_plink.sh \                               # main script
../../../toy_data/genotype \                            # genotype data in plink format
../../../toy_data/pheno_continue.txt  \                 # pheno file
SEX,AGE \                                               # covariate list
PHENO \                                                 # target variable
quantitative                                            # type target variable



```

### 🔧 GMMAT Pipeline

**Run with:**

```bash

1. Get inside the GMMAT folder.
2. Run the command:
      ./run_pipeline_gmmat.sh ./ ../../../toy_data/input_kinship ../../../toy_data/input_dosage ../../../toy_data/pheno_binary.txt "PHENO ~ AGE + SEX" quantitative

# explanation

./run_pipeline_gmmat.sh  \                      # the main script
./ \                                            # location of the working directory
../../../toy_data/input_kinship \               # location of the genotype data for kinship
../../../toy_data/input_dosage \                # location of the genotype (dosage)
../../../toy_data/pheno_binary.txt \            # location of pheno file
"PHENO ~ AGE + SEX" \                           # model selected
quantitative                                    # type of pheno variable (quantitative or binary)


```

### 🔧 SAIGE Pipeline

**Run with:**

```bash

1. Get inside the SAIGE folder.
2. Run the command:
      cp /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/SAIGE/test_our_pipeline/SINGLE_MARKER/Saige_1.3.0.sif ./
3. Run the command
     ./run_pipeline_saige.sh ../../../toy_data/input_kinship ../../../toy_data/input_dosage ../../../toy_data/pheno_binary.txt  AGE,SEX, SEX PHENO quantitative

# Explanation

./run_pipeline_saige.sh \                       # main script
../../../toy_data/input_kinship \               # kinship input
../../../toy_data/input_dosage \                # dosage input
../../../toy_data/pheno_binary.txt  \           # pheno data
AGE,SEX, \                                      # list of covariates
SEX \                                           # binary covariates
PHENO \                                         # target variable
quantitative                                    # type of variable


```

- Note: Saige require the image: Saige_1.3.0.sif inside the SAIGE folder. Please, just copy the saige image from my folder to your working directory: 

  cp /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/SAIGE/test_our_pipeline/SINGLE_MARKER/Saige_1.3.0.sif ./

##  📊 Outputs
Each pipeline will generate a folder with the following information:

✅ GWAS summary statistics (sum_stat.txt)

📈 Manhattan plot (Rect_Manhtn.PVAL.jpg)

📉 QQ plot (qqplot_with_lambda.jpg)



