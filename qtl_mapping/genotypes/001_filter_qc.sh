#!/bin/bash
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}

module purge
module load bcftools
module load tabix

ROOTDIR='/gcglm/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf'
OUTDIR='/gcglm/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf_AF5'

mkdir -p /gcglm/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf_AF5/stat_reports

echo "Filtering variants >5%"

#bcftools view -i 'RAF > 0.05' "$ROOTDIR/chr${i}_ashkenazi.407.vcf.gz" | bcftools view -i 'RAF < 0.99' -Oz -o "$OUTDIR/chr${i}_ashkenazi.407.AF5.vcf.gz" --threads 1

#bcftools index --tbi "$OUTDIR/chr${i}_ashkenazi.407.AF5.vcf.gz"
echo "Filtering finished"

#bcftools stats "$OUTDIR/chr${i}_ashkenazi.407.AF5.vcf.gz" > $OUTDIR/stat_reports/chr${i}_AF5.stat
#bcftools view "$OUTDIR/chr${i}_ashkenazi.407.AF5.vcf.gz" | grep -v "^#" | wc -l
echo "Filtering >99% call rate"

#bcftools view -i 'F_MISSING < 0.05' "$OUTDIR/chr${i}_ashkenazi.407.AF5.vcf.gz" -Oz -o "$OUTDIR/chr${i}_ashkenazi.407.AF5.QC5.vcf.gz" --threads 1

#bcftools index --tbi "$OUTDIR/chr${i}_ashkenazi.407.AF5.QC5.vcf.gz"
echo "Quality control finished"

bcftools stats "$OUTDIR/chr${i}_ashkenazi.407.AF5.QC5.vcf.gz" > $OUTDIR/stat_reports/chr${i}_QC5.stat
bcftools view "$OUTDIR/chr${i}_ashkenazi.407.AF5.QC5.vcf.gz" | grep -v "^#" | wc -l
