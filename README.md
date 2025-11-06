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
  - [Singularity 4.2.1] (https://docs.sylabs.io/guides/latest/user-guide.pdf)
  - [R 4.2.2](https://cran.r-project.org/)

- **R Packages**
  - [GMMAT](https://cran.r-project.org/web/packages/GMMAT/)
  - [ggplot2](https://ggplot2.tidyverse.org/)
  - [qqman](https://cran.r-project.org/web/packages/qqman/)

---

## 5. Installation

```bash

git clone https://github.com/bciezah1/SAGA.git

cd SAGA/single_marker/tools/bin
chmod +x *

cd SAGA/single_marker/tools/PLINK
chmod +x *

cd SAGA/single_marker/tools/GMMAT
chmod +x *

cd SAGA/single_marker/tools/SAIGE
chmod +x *

```
---

## 6. Repository Structure

```
SAGA/
├── single_marker/ # Main pipeline scripts and tools
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

        ./run_pipeline_plink.sh full_path_to_geno/geno  full_path_to_pheno/pheno_binary.txt  COV1,COV2,PC1,PC2,PC3  PHENO  binary  myoutputs


        # Explanation

        ./run_pipeline_plink.sh \                               # main script
        full_path_to_geno/geno  \                               # genotype data in plink format
        full_path_to_pheno/pheno_binary.txt  \                  # pheno file
        COV1,COV2,PC1,PC2,PC3 \                                 # covariate list (up to 5)
        PHENO \                                                 # target variable
        binary                                                  # type target variable
        myoutputs                                               # name of your folder where results will be located

    


```

### 🔧 GMMAT Pipeline

**Run with:**

```bash

1. Get inside the GMMAT folder.
2. Run the command:

      ./run_pipeline_gmmat.sh  full_path_to_geno/geno  full_path_to_pheno/pheno_binary.txt  "PHENO ~ COV1 + COV2 + COV3 + PC1 + PC2 + PC3"  binary  myoutput
      
      # explanation
      ./run_pipeline_gmmat.sh  \                      # the main script
      full_path_to_geno/geno  \                       # genotype data in plink format
      full_path_to_pheno/pheno_binary.txt  \          # pheno file
      "PHENO ~ COV1 + COV2" \                         # model selected
      binary \                                        # type of pheno variable (quantitative or binary)
      myoutput                                        # location of my output


```

### 🔧 SAIGE Pipeline

**Run with:**

```bash

1. Download the Saige image using this link: https://github.com/bciezah1/SAGA/releases/download/v1.0/Saige_1.3.0.sif to get the saige image, and put the image in the SAIGE folder

2. Get inside the SAIGE folder.

3. Run the command

      ./run_pipeline_saige.sh  full_path_to_geno/geno  full_path_to_pheno/pheno_binary.txt   COV1,COV2,PC1,PC2,PC3  COV1  PHENO  binary  myoutput
      
      # Explanation
      
      ./run_pipeline_saige.sh \                       # main script
      full_path_to_geno/geno  \                       # genotype data in plink format
      full_path_to_pheno/pheno_binary.txt  \          # pheno file
      COV1,COV2,PC1,PC2,PC3,PC4,PC5 \                 # list of covariates
      COV1 \                                          # binary covariates
      PHENO \                                         # target variable
      binary \                                        # type of variable
      myoutput                                        # working directory


```
## 8. Input Formats

SAGA only need *two* inputs, (1) genetic data in plink format (`.bed`, `.bim`, `.fam`), and a (2) phenotype file.  

---

### **8.1. Genetic data**  
In plinnk format
  
**Example:**

```bash
    geno.bed
    geno.bim
    geno.fam
```

  > 💡 You can input your raw genetic data, and SAGA will take care of the QC steps.

---

### **8.2. Phenotype File**  
A tab-delimited file with sample IDs, phenotype values, and optional covariates.

**Required columns:**

FID  IID  PHENO

**Recommended:**

FID  IID  COV1  COV2  PHENO

**Example:**

```bash

FID	    IID	    COV1	COV2	PHENO
FAM001	IND001	0	    84	    1
FAM002	IND002	0	    85	    0
FAM003	IND003	1	    72	    1


...

```
### **📌 How They Fit in the Workflow**

![Pipeline Diagram](https://github.com/bciezah1/SAGA/blob/main/images/diagram.jpg)

> ⚠️ **Warning:** The number of participants and their order **MUST match exactly** between the PLINK files and the phenotype file.

> ⚠️ **Important:** Refer to the `toy_data/` folder to verify correct formatting and file naming.

>  ⚠️ **Important:** Open permit to executables open permits for **ALL** files on each folder

       Tip: To open permits, just go inside the folder mentioned, and run the next command: chmod +x *




---

## 9.  📊 Output Examples
Each pipeline will generate a folder with the following information:

✅ GWAS summary statistics (sum_stat.txt)

```bash
SNP                     CHR     POS     Allele1 Allele2 AF      BETA            SE               PVAL    
chr1:2917460:G:A        1       2917460 G       A       0.364   0.19181         0.111096        8.425454E-02    
chr1:3064229:T:C        1       3064229 T       C       0.09    0.154274        0.184408        4.028223E-01    
chr1:3086968:G:A        1       3086968 G       A       0.169   -0.145269       0.144045        3.132143E-01    
chr1:3173790:G:A        1       3173790 G       A       0.0425  -0.069691       0.269407        7.958797E-01    
chr1:3447923:C:T        1       3447923 C       T       0.063   -0.171615       0.222336        4.401907E-01    
chr1:3906831:C:T        1       3906831 C       T       0.2045  -0.0609915      0.133792        6.484836E-01    
chr1:3949011:T:C        1       3949011 T       C       0.222   0.121489        0.134684        3.670415E-01    
chr1:4109988:C:T        1       4109988 C       T       0.2065  -0.151731       0.130497        2.449444E-01    
chr1:4474909:C:T        1       4474909 C       T       0.184   0.29204         0.13815         3.452163E-02    

...

```

📈 Manhattan plot 

![Pipeline Diagram](https://github.com/bciezah1/SAGA/blob/main/images/rectangular_manhattan_plot.jpg)

📉 QQ plot 

![Pipeline Diagram](https://github.com/bciezah1/SAGA/blob/main/images/Q_Q_plot.jpg)

📈 SNP density plot

![Pipeline Diagram](https://github.com/bciezah1/SAGA/blob/main/images/SNP_density_plot.jpg)

📉 Circular Manhattan plot

![Pipeline Diagram](https://github.com/bciezah1/SAGA/blob/main/images/circular_manhattan_plot.jpg)



