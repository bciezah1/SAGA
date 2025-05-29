# EeasyGWAS 

## 🧾 Summary

EasyGWAS is a collection of pipelines to perform single marker GWAS using GMMAT and SAIGE. This tool has been developed with the goal of allowing no programmer scientist to perform GWAS on their genotyping/imputated data. 

## 🧾 Required Inputs

EasyGWAS require 3 main inputs (listed below). To reduce issues, you should rename your corresponding files as they are label on our toy_data folder.

        - High quality imputation data (plink format) to calculate kinship (R2=0.99 & MAF=0.05)
        - Dosage data (plink format) to perform association (Suggested R2=0.80 & MAF=0.01)
        - Pheno file (Should include FID, IID, ADRD, SEX, AGE, and APOEe4).

        Note: To check format and name of the files, please, review the format/label/etc on the toy_data folder.

## 🧾 Pipelines

Two main pipelines have been implemented (GMMAT & SAIGE) so far. We recommend to first run the pipelines (as suggested bellow), and then replace your files by the corresponding files on the toy_data and re run the same command. 
        
        ## GMMAT

                ## How to Run

                - bash submit_all.sh ./ ./../toy_data/ ../toy_data/ model1

        ## SAIGE

                ## How to Run

                - sbatch saige.slurm ../toy_data/ model1
## 🧾 Outputs:

         - Summary Statistics
         - Manhattan plot
         - QQplot


