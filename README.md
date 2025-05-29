# EasyGWAS

## 🧬 Summary

**EasyGWAS** is a collection of streamlined pipelines for performing single-marker GWAS using **GMMAT** and **SAIGE**. It is designed for scientists with limited programming experience who want to analyze genotyped or imputed genetic data with minimal setup.

---

## 📥 Required Inputs

EasyGWAS requires **three main input files**. To minimize errors, we recommend renaming your files to match the naming conventions shown in the `toy_data/` folder.

### Input Files:

- **Kinship input (PLINK format):** High-quality imputed genotype data for kinship estimation.  
  _Recommended: `R² ≥ 0.99`, `MAF ≥ 0.05`_

- **Association input (PLINK format):** Dosage or genotype data for association testing.  
  _Recommended: `R² ≥ 0.80`, `MAF ≥ 0.01`_

- **Phenotype file:** Must include the following columns:  
  `FID`, `IID`, `ADRD`, `SEX`, `AGE`, `APOEe4`

> 📎 _Refer to the `toy_data/` folder to verify correct formatting and file naming._

---

## 🚀 Pipelines

EasyGWAS includes two pipelines: one using **GMMAT** and one using **SAIGE**. We recommend running the provided toy data first to verify that everything is working before applying the pipeline to your own dataset.

### 🔧 GMMAT Pipeline

**Run with:**

```bash
bash submit_all.sh ./ ../toy_data/ ../toy_data/ model1

```

### 🔧 SAIGE Pipeline

**Run with:**

```bash
sbatch saige.slurm ../toy_data/ model1
```

- Note: Saige require the image: Saige_1.3.0.sif inside the SAIGE folder. Please, just copy the saige image from my folder to your working directory: 

  cp /mnt/vast/hpc/gtosto_lab/GT_ADMIX/Basilio_08_19_2022/GWAS/SAIGE/test_our_pipeline/SINGLE_MARKER/Saige_1.3.0.sif ./

##  📊 Outputs
Each pipeline will generate the following:

✅ GWAS summary statistics

📈 Manhattan plot

📉 QQ plot



