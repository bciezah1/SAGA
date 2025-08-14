# SAGA

## 1. Summary

**SAGA** (Simplified Association Genomewide Analyses) is a collection of streamlined pipelines for performing single-marker GWAS using **PLINK**, **GMMAT** and **SAIGE**. It is designed for scientists with limited programming experience who want to analyze genotyped or imputed genetic data with minimal setup.

---

## 2. Background

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

## 3. Features

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

## 4. Dependencies

- **Operating System**
  - UNIX/Linux environment

- **Software**
  - [R 4.2.2](https://cran.r-project.org/)

- **R Packages**
  - [GMMAT](https://cran.r-project.org/web/packages/GMMAT/)
  - [ggplot2](https://ggplot2.tidyverse.org/)
  - [qqman](https://cran.r-project.org/web/packages/qqman/)

---

## 5. Installation

```bash

git clone https://github.com/bciezah1/SAGA.git
cd SAGA

```
---

## 6. Repository Structure

```
SAGA/
├── SINGLE_MARKER/ # Main pipeline scripts and tools
│ ├── tools/ # Tool-specific scripts and binaries
│ │ ├── bin/ # Binaries for GWAS tools
│ │ ├── GMMAT/ # GMMAT workflow scripts
│ │ ├── PLINK/ # PLINK workflow scripts
│ │ └── SAIGE/ # SAIGE workflow scripts and container (https://github.com/bciezah1/SAGA/releases/download/v1.0/Saige_1.3.0.sif)
└── toy_data/ # Example genotype and phenotype data

```
---

##  7. How to Run it

SAGA includes three pipelines: The first one using **PLINK**, the second one using **GMMAT**, and third one using **SAIGE**. We recommend running the provided toy data first to verify that everything is working properly before applying the pipeline to your own dataset.

### 🔧 PLINK Pipeline

**Run with:**

```bash

1. Get inside the PLINK folder.
2. Run the command:
      ./run_pipeline_plink.sh ../../../toy_data/genotype ../../../toy_data/pheno_continue.txt  COV1,COV2,COV3 PHENO quantitative


# Explanation

./run_pipeline_plink.sh \                               # main script
../../../toy_data/genotype \                            # genotype data in plink format
../../../toy_data/pheno_continue.txt  \                 # pheno file
COV1,COV2,COV3 \                                        # covariate list
PHENO \                                                 # target variable
quantitative                                            # type target variable


```

### 🔧 GMMAT Pipeline

**Run with:**

```bash

1. Get inside the GMMAT folder.
2. Run the command:
      ./run_pipeline_gmmat.sh ./ ../../../toy_data/input_kinship ../../../toy_data/input_dosage ../../../toy_data/pheno_continue.txt "PHENO ~ COV1 + COV2 + COV3" quantitative

# explanation

./run_pipeline_gmmat.sh  \                      # the main script
./ \                                            # location of the working directory
../../../toy_data/input_kinship \               # location of the genotype data for kinship
../../../toy_data/input_dosage \                # location of the genotype (dosage)
../../../toy_data/pheno_binary.txt \            # location of pheno file
"PHENO ~ COV1 + COV2 + COV3" \                  # model selected
quantitative                                    # type of pheno variable (quantitative or binary)


```

### 🔧 SAIGE Pipeline

**Run with:**

```bash

1. Download the Saige image from https://github.com/bciezah1/SAGA/releases/download/v1.0/Saige_1.3.0.sif and put the image in the SAIGE folder

2. Get inside the SAIGE folder.

2. Run the command
     ./run_pipeline_saige.sh ../../../toy_data/input_kinship ../../../toy_data/input_dosage ../../../toy_data/pheno_continue.txt  COV1,COV2,COV3 COV1 PHENO quantitative

# Explanation

./run_pipeline_saige.sh \                       # main script
../../../toy_data/input_kinship \               # kinship input
../../../toy_data/input_dosage \                # dosage input
../../../toy_data/pheno_binary.txt  \           # pheno data
COV1,COV2,COV3 \                                # list all covariates (quantitative + binary)
COV1 \                                          # only binary covariates
PHENO \                                         # target variable
quantitative                                    # type of variable


```
## Input formats

- **Kinship input (PLINK format):** High-quality imputed genotype data for kinship estimation.  

  expected label file:
  
      input_kinship.bed
      input_kinship.bim
      input_kinship.fam

- **Dosage input (PLINK format):** Dosage or genotype data.  

  expected label files:

      input_dosage.bed
      input_dosage.bim
      input_dosage.fam

- **Phenotype file:** 
  
  Must include the following columns:  
  
  `FID`, `IID`, `PHENO`

  recommended:

  `FID`, `IID`, `PHENO`, `COV1`, `COV2`, `COV3`, `COV4`, `COV5` 
  
  example pheno file format:

      FID     IID     COV1    COV2    COV3    COV4    COV5    PHENO
      FAM001  IND001  0       84      9.24046 5.93909 3.06394 4.87729133902262
      FAM002  IND002  0       85      5.78941 7.40133 7.86926 9.76982251051672
      FAM003  IND003  1       72      4.3637  3.32195 7.7888  2.35685104797102
      FAM004  IND004  0       69      1.00887 7.85084 8.35159 4.90705898152053
      FAM005  IND005  1       77      7.61209 4.96077 4.26298 5.43028469635999


      ...


> ⚠️ **Warning:** The number of participants and their order **MUST match exactly** between the PLINK files and the phenotype file.

> ⚠️ **Important:** Refer to the `toy_data/` folder to verify correct formatting and file naming.

>  ⚠️ **Important:** Open permit to executables
       Open permits for *ALL* files on each folder
       Tip: To open permits, just go inside the folder mentioned, and run the next command: chmod +x *


---

##  📊 Outputs
Each pipeline will generate a folder with the following information:

✅ GWAS summary statistics (sum_stat.txt)

📈 Manhattan plot (Rect_Manhtn.PVAL.jpg)

📉 QQ plot (qqplot_with_lambda.jpg)



