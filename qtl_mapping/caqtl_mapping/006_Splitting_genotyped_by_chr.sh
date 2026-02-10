#!/bin/bash
#SBATCH --job-name=generatePLINKs                # Job name
#SBATCH --partition=pe2                     # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail
#SBATCH --mem=20G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=24:00:00                       # Time limit 4 hours
#SBATCH --output=logs/005_splitting_genotyped_by_chr_%A.log               # Standard output and error log
#SBATCH --array=1-22

module load miniconda3/3.22.0
source activate genomic_tools #preinstall jupyter in the conda environment in which your job will be conducted

plink_file="/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_101424_cpm_tmm_maf5_fdr5_50kb/results/004_genotypes/plink/CD4_all_chr_ashkenazi.362.AF1.QC.BA.king2.hwe.annot"  #this big plink file is the result of concatenating all preprocessed chr_vcfs from Sam's genotyping analysis
#sam's analysis can be found at /gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf
#this file has gone through various phases of QC for whioch samples were removed from the original dataset, based on outliers and common samples between WGS and GEX 


OUTDIR="/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_101424_cpm_tmm_maf5_fdr5_50kb/results/004_genotypes/plink"

mkdir -p $OUTDIR/per_chr_plink_files

plink2 --bfile $plink_file   \
        --chr ${SLURM_ARRAY_TASK_ID} \
        --make-bed --out $OUTDIR/per_chr_plink_files/chr${SLURM_ARRAY_TASK_ID}_ashkenazi.362.AF1.QC.BA.king2.hwe.annot