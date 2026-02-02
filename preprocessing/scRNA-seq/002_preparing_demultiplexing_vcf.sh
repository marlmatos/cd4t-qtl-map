#!/bin/bash -l
#SBATCH --job-name=bcftools     # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=marliette.rodriguezmatos@einsteinmed.edu    
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=32G                     # Job memory request
#SBATCH --output=002_preparing_demultiplexing_vcf-%j.out

module load bcftools/1.19

data_dir='/gcgl/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf'
dir='/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf'

# subset all chromosomes for just sample 30193 ( a random sample, what we want is just that column)

echo "Iterating over chromosome files"
for i in {1..22}; do
    chr_file="${data_dir}/chr${i}_ashkenazi.407.AF1.QC.BA.vcf.gz"
    # Get the filename without extension
    filename=$(basename "$chr_file" .vcf.gz)
    # Uncompress VCFs and output to temporary VCF
    bcftools view "$chr_file" -Ou -o "temp_${filename}.vcf"
done

ls temp_*.vcf > vcf_path_list.txt

echo "Concatenating all chromosome files"
bcftools concat -f vcf_path_list.txt -Oz -o $dir/CD4_allsamples_common_maf1.vcf.gz 

echo "Removing temporary files"
rm temp_*