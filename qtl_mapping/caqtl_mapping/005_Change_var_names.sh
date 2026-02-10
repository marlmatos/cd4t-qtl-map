#!/bin/bash
#SBATCH --job-name=plink        # Job name
#SBATCH -p pe2
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G                     # Job memory reques
#SBATCH --mail-type=END,FAIL             # send email when job begins
#SBATCH --mail-user=mmatos@nygenome.org
#SBATCH --time=24:00:00             # Wall time limit (days-hrs:min:sec)
#SBATCH --output=logs/005_change_var_names.out    # Path to the standard output and error files relative to the working

OUTDIR="/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/004_genotypes/plink"

module load miniconda3/3.22.0
source activate genomic_tools #preinstall jupyter in the conda environment in which your job will be conducted


plink2 --bfile $OUTDIR/CD4_all_chr_ashkenazi.362.AF1.QC.BA.king2.hwe \
--memory 12000 \
--set-all-var-ids @:#[b38]\$r,\$a \
--new-id-max-allele-len 198 \
--make-bed --out $OUTDIR/CD4_all_chr_ashkenazi.362.AF1.QC.BA.king2.hwe.annot 





