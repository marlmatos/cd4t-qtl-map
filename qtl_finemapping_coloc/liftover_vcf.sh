#!/bin/bash
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10

module load bcftools

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz{,.tbi} bcftools +liftover --no-version -Ou 
ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz -- \
  -s $HOME/GRCh37/human_g1k_v37.fasta \
  -f $HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c $HOME/GRCh38/hg19ToHg38.over.chain.gz \
  --reject ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.reject.bcf \
  --reject-type b \
  --write-src | \ bcftools sort -Ob | tee ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.hg38.bcf | \
bcftools index --force --output ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.hg38.bcf.csi
