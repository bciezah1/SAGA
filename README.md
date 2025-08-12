# SAGA

- ## Summary

**SAGA** is a collection of streamlined pipelines for performing single-marker GWAS using **PLINK**, **GMMAT** and **SAIGE**. It is designed for scientists with limited programming experience who want to analyze genotyped or imputed genetic data with minimal setup.

---

- ## Background

A **Genome-Wide Association Study (GWAS)** is a statistical method used to scan the entire genome for genetic variants—most often single nucleotide polymorphisms (SNPs)—that are associated with specific traits or diseases.  
By comparing genetic data from many individuals, GWAS can help identify variants linked to:
- Common diseases (e.g., diabetes, Alzheimer’s)
- Physical traits (e.g., height, cholesterol levels)
- Drug responses and adverse reactions

#### Why GWAS matters
GWAS has transformed our understanding of the genetic architecture of complex traits, providing clues for:
- Risk prediction
- Disease prevention
- Development of targeted therapies

#### The problem
Despite their value, GWAS can be **technically challenging** for newcomers because:
- Data needs careful quality control and preprocessing.
- Multiple specialized tools are required.
- Analyses often require complex scripting and computational environments.

#### The solution: SAGA
**SAGA** simplifies this process into a **single, automated pipeline**.  

Users provide:

- Genotype data (in PLINK format)
- Phenotype file (traits of interest)
  
SAGA handles:

- Quality control  
- Population structure adjustment  
- Choice of optimal GWAS tool (PLINK, GMMAT, or SAIGE)  
- Result visualization (Manhattan + QQ plots)  

This lowers the entry barrier for clinicians, researchers, and students without advanced computational skills.

---

- ## Features

- **Three integrated GWAS backends:**
  - **PLINK** — fast linear/logistic regression for unrelated individuals.
  - **GMMAT** — mixed models for family/related samples.
  - **SAIGE** — scalable logistic mixed models for large and unbalanced case-control datasets.
- **Automated preprocessing:**
  - Sample/variant QC  
  - PCA for population structure adjustment  
  - Kinship matrix generation
- **Standardized outputs:**
  - Summary statistics
  - High-resolution plots (`.png`, `.pdf`)
- **Fully bash-based:** no SLURM or HPC dependency, works on any UNIX/Linux environment.

---

- ## Dependencies

- **Operating System**
  - UNIX/Linux environment

- **Software**
  - [R 4.2.2](https://cran.r-project.org/)

- **R Packages**
  - [GMMAT](https://cran.r-project.org/web/packages/GMMAT/)
  - [ggplot2](https://ggplot2.tidyverse.org/)
  - [qqman](https://cran.r-project.org/web/packages/qqman/)

---

- ## Installation

```bash

git clone https://github.com/bciezah1/SAGA.git
cd SAGA

```

- ## Repository Structure

```
SAGA/
├── SINGLE_MARKER/ # Main pipeline scripts and tools
│ ├── tools/ # Tool-specific scripts and binaries
│ │ ├── bin/ # Binaries for GWAS tools
│ │ ├── GMMAT/ # GMMAT workflow scripts
│ │ ├── PLINK/ # PLINK workflow scripts
│ │ └── SAIGE/ # SAIGE workflow scripts and container
└── toy_data/ # Example genotype and phenotype data

```

## 📥 Required Inputs

SAGE requires **three main input files**. To minimize errors, we recommend renaming your files to match the naming conventions shown in the `toy_data/` folder.

### Input Files:

- **Kinship input (PLINK format):** High-quality imputed genotype data for kinship estimation.  
  _Recommended: `R² ≥ 0.99`, `MAF ≥ 0.05`_

  expected label file:
  
      input_kinship.bed
      input_kinship.bim
      input_kinship.fam

- **Dosage input (PLINK format):** Dosage or genotype data for association testing.  
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

2. Run the command
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


##  📊 Outputs
Each pipeline will generate a folder with the following information:

✅ GWAS summary statistics (sum_stat.txt)

📈 Manhattan plot (Rect_Manhtn.PVAL.jpg)

📉 QQ plot (qqplot_with_lambda.jpg)



