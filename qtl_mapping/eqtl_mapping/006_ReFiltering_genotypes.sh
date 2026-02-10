#!/bin/bash
#SBATCH --job-name=plink        # Job name
#SBATCH -p cpu
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --time=24:00:00             # Wall time limit (days-hrs:min:sec)
#SBATCH --output=006_filt_format.out    # Path to the standard output and error files relative to the working

###################################################
# # **Title:** Filtering geneotypes to match samples that pass all filters (samples that failed in scRNAseq but passed in WGS)
# # **Author:** Marliette Matos
# # **Date:** 08/16/2024
# # **Description:** Using genotype files QCed by Sam Gathan and filtered for a MAF>5%

module load miniconda3/3.22.0
source activate genomic_tools #preinstall jupyter in the conda environment in which your job will be conducted


#Enviroment variables
OUTDIR="/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5"
COMMON_SAMPLES="/gpfs/commons/home/mmatos/cd4_QTL_analysis/02_Gene_expression/analysis/001_preparing_seurat_obj_QTL/003_common_samples_gfgex_wgs.in.tsv"
#common samples between WGS and Gene Expression

# Keeping only common samples  
plink2 --bfile $OUTDIR/CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe --memory 12000 --keep $COMMON_SAMPLES --make-bed --out $OUTDIR/CD4_all_chr_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq

# Before piping the resulting into tensorqTL, the variants should be renamed by chr/post/ref/alt -> 003_change_var_names.sh
plink2 --bfile $OUTDIR/CD4_all_chr_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq \
--memory 12000 \
--set-all-var-ids @:#[b38]\$r,\$a \
--new-id-max-allele-len 198 \
--make-bed --out $OUTDIR/CD4_all_chr_ashkenazi.363.AF1.QC.BA.king2.hwe.scrnaseq.annot 

# Splittting the pre-pruned by chromosomes for cis-eQTL calling -> 007_split_plinkfiles_by_chr.sh