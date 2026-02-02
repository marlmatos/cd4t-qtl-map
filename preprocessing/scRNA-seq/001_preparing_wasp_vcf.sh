#!/bin/bash -l
#SBATCH --job-name=bcftools     # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=marliette.rodriguezmatos@einsteinmed.edu    
#SBATCH -p pe2 
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=4                   
#SBATCH --mem=32G                     # Job memory request
#SBATCH --output=001_preparing_wasp_vcf-%j.out

module load bcftools/1.19

data_dir='/gcgl/mmatos/cd4_aging_project/data/LowPass_WGS/03.27.24_SamGhatan_QC/ashkenazi_vcf'
dir='/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf'

# subset all chromosomes for just sample 30193 ( a random sample, what we want is just that column)

echo "Iterating over chromosome files"
for i in {1..22}; do
    chr_file="${data_dir}/chr${i}_ashkenazi.407.AF1.QC.BA.vcf.gz"
    # Get the filename without extension
    filename=$(basename "$chr_file" .vcf.gz)
    # Extract sample 10013 and output to temporary VCF
    bcftools view -v snps -s 10013 "$chr_file" -Ou -o "temp_${filename}.vcf"
done

ls temp_*.vcf > vcf.list
echo "Concatenating over chromosome files"

bcftools concat -f vcf.list -Ou -o $dir/temp_1.vcf 


echo "Clearing the info,ID and Format field because the tags"

bcftools annotate --remove INFO,FORMAT/RC,FORMAT/AC,FORMAT/DP,FORMAT/DS,FORMAT/GP -Ou -o $dir/temp_2.vcf $dir/temp_1.vcf

echo "Changing all genotype instances to just (0/1)"
sed 's/\(0\/1\|0\/0\|1\/0\|1\/1\)/0\/1/g' $dir/temp_2.vcf > $dir/temp_3.vcf

sed 's/,length=[0-9]\+//g' $dir/temp_3.vcf > $dir/temp_4.vcf

echo "Removing duplicate variants"
bcftools norm --rm-dup all $dir/temp_4.vcf -Ou -o $dir/CD4p_WGS_pass_only.snps.maf1.WASP.vcf 


bcftools stats $dir/CD4p_WGS_pass_only.snps.maf1.WASP.vcf > $dir/CD4p_WGS_pass_only.snps.maf1.WASP.stats

echo "Indexing final vcf"
bcftools index -t $dir/CD4p_WGS_pass_only.snps.maf1.WASP.vcf 

echo "Removing temporary files"
rm temp_*
