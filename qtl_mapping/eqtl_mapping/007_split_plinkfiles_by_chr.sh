#!/bin/bash
#SBATCH --job-name=generatePLINKs                # Job name
#SBATCH --partition=cpu                     # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail
#SBATCH --mem=20G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=24:00:00                       # Time limit 4 hours
#SBATCH --output=logs/006_splitting_genotyped_by_chr_%A_%a.log               # Standard output and error log
#SBATCH --error=logs/006_splitting_genotyped_by_chr_%A_%a.err
#SBATCH --array=1-22

###################################################
# # **Title:** Splint splink files per chromosomes
# # **Author:** Marliette Matos
# # **Date:** 08/16/2024

module load miniconda3/3.22.0
source activate genomic_tools

plink_file="/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq.annot"  #this big plink file is the result of concatenating all preprocessed chr_vcfs from Sam's genotyping analysis
#sam's analysis can be found at /gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf
#this file has gone through various phases of QC for whioch samples were removed from the original dataset, based on outliers and common samples between WGS and GEX 

OUTDIR="/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/"

mkdir -p $OUTDIR/per_chr_plink_files_scRNAseq

plink2 --bfile $plink_file   \
        --chr ${SLURM_ARRAY_TASK_ID} \
        --make-bed --out $OUTDIR/per_chr_plink_files_scRNAseq/chr${SLURM_ARRAY_TASK_ID}_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq.annot