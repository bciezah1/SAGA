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
        R/4.2.2
>         
        PLINK/1.9.10


---

## 🚀 Pipelines

EasyGWAS includes two pipelines: one using **GMMAT** and one using **SAIGE**. We recommend running the provided toy data first to verify that everything is working before applying the pipeline to your own dataset.

### 🔧 GMMAT Pipeline

**Run with:**

```bash

1. Get inside the GMMAT folder.
2. Run the command:
      bash submit_all.sh ./ ../toy_data/ ../toy_data/ model1

```

### 🔧 SAIGE Pipeline

**Run with:**

```bash

1. Get inside the SAIGE folder.
2. Run the command:
    cp /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/SAIGE/test_our_pipeline/SINGLE_MARKER/Saige_1.3.0.sif ./
3. Run the command
    sbatch saige.slurm ../toy_data/ model1

```

- Note: Saige require the image: Saige_1.3.0.sif inside the SAIGE folder. Please, just copy the saige image from my folder to your working directory: 

  cp /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/SAIGE/test_our_pipeline/SINGLE_MARKER/Saige_1.3.0.sif ./

##  📊 Outputs
Each pipeline will generate the following:

✅ GWAS summary statistics (sum_stat.txt)

📈 Manhattan plot (Rect_Manhtn.PVAL.jpg)

📉 QQ plot (qqplot_with_lambda.jpg)



