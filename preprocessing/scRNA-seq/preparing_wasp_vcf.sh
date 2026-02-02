#!/bin/bash -l
#SBATCH --job-name=bcftools     # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=marliette.rodriguezmatos@einsteinmed.edu    
#SBATCH -p normal 
#SBATCH -t 2-00:00:00
#SBATCH --nodes=2                    # Run on a single CPU
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --output=bcftools-%j.out

conda activate genomic_tools

fasta=/ggum/genome/Gencove_WGS/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

#plink2 --vcf merged_CD4p_WGS_pass_only.vcf.gz \
 #   --const-fid 0 \
 #       --snps-only 'just-acgt' \
 #           --memory 120000 \
 #               --make-bed \
 #                   --out CD4p_WGS_pass_only.snps

#plink2 --bfile CD4p_WGS_pass_only.snps \
 #   --maf 0.05 \
 #   --make-bed \
 #   --out CD4p_WGS_pass_only.snps.maf05 

#plink2 --bfile CD4p_WGS_pass_only.snps.maf05 \
 #   --export vcf \
 #   --out CD4p_WGS_pass_only.snps.maf05

dir=/ggum/aging_project/ATAC-seq_analysis/scripts/nf-core_atac/try_11072023_star_allele_specfic/var_vcf
plink2 --bfile CD4p_WGS_pass_only.snps.maf05 \
    --fa $fasta --ref-from-fa \
    --export vcf \
    --out $dir/CD4p_WGS_pass_only.snps.maf05

conda deactivate

conda activate bcf_env

bcftools view -s 30193  $dir/CD4p_WGS_pass_only.snps.maf05.vcf > $dir/temp_CD4p_WGS_pass_only.snps.maf05.v1.vcf

#clearing the info field because the tags
bcftools annotate --remove INFO -o $dir/temp_CD4p_WGS_pass_only.snps.maf05.v2.vcf $dir/temp_CD4p_WGS_pass_only.snps.maf05.v1.vcf


#add "chr" to all chromosomes
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $dir/temp_CD4p_WGS_pass_only.snps.maf05.v2.vcf > $dir/temp_CD4p_WGS_pass_only.snps.maf05.v2chr.vcf

#change all genotype instances to just 
sed 's/\(0\/1\|0\/0\|1\/0\|1\/1\)/0\/1/g' $dir/temp_CD4p_WGS_pass_only.snps.maf05.v2chr.vcf > $dir/temp_CD4p_WGS_pass_only.snps.maf05.v3chr.vcf

#adds "chr" to the contig name in the header
sed -i 's/##contig=<ID=/##contig=<ID=chr/g' $dir/temp_CD4p_WGS_pass_only.snps.maf05.v3chr.vcf

sed 's/,length=[0-9]\+//g' $dir/temp_CD4p_WGS_pass_only.snps.maf05.v3chr.vcf > $dir/CD4p_WGS_pass_only.snps.maf05.v4chr.vcf

bcftools norm --rm-dup all $dir/CD4p_WGS_pass_only.snps.maf05.v4chr.vcf -Oz -o $dir/CD4p_WGS_pass_only.snps.maf05.v5chr.vcf.gz 

#bcftools view $dir/CD4p_WGS_pass_only.snps.maf05.v4chr.vcf -Oz -o $dir/CD4p_WGS_pass_only.snps.maf05.v5chr.vcf.gz 
bcftools index -t $dir/CD4p_WGS_pass_only.snps.maf05.v5chr.vcf.gz

rm temp_*