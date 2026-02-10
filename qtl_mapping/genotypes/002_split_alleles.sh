#!/bin/bash
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22
#SBATCH --output=logs/002_split_alleles-%j.out

i=${SLURM_ARRAY_TASK_ID}

module purge
module load bcftools
module load tabix

DATA='/gcglmcd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf'
OUTDIR="/gcglmcd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf_AF5"

bcftools norm -m-any "$DATA/chr${i}_ashkenazi.407.vcf.gz" | bcftools view -i 'RAF > 0.05' | bcftools view -i 'RAF < 0.99' | bcftools view -i 'F_MISSING < 0.05' -Oz -o "$OUTDIR/chr${i}_ashkenazi.407.AF5.Q5.BA.vcf.gz"

bcftools stats "$OUTDIR/chr${i}_ashkenazi.407.AF5.Q5.BA.vcf.gz" > $OUTDIR/stat_reports/chr${i}_QC5.BA.stat
