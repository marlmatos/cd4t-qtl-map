#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# 01. Define regions containing GWS SNPs
# Provide eQTL sumstats path
# eQTL name

#2
numbers=($(seq 1 22))

for i in "${numbers[@]}";do

echo "Running chr${i}"

module purge
module load R/4.4.1

NAME="CD4T_chromatin"
DIR=$(pwd)
EQTL="/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr${i}.csv"
INDEP="/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/006_caQTLs/filtered_qsmooth_norm_cpm_1mb/cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr${i}.csv"

Rscript 01.Define_finemap_regions.caQTL.R \
	"${INDEP}" \
	"${NAME}" \
	"${DIR}"
	

# Calculate LD for regions containing significant SNPs
sbatch --wait 02.Generate_LD.sh "regions/${NAME}.chr${i}.regions.txt" \
			 "/gchm/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/per_chr_plink_files/chr${i}_ashkenazi.367.AF1.QC.BA.king2.hwe.annot" \
			 "Ashkenazi"

# Wait for job to finish
wait

echo "Running finemapping"

# Finemap regions
sbatch 03.submit_caQTL_finemap.sh \
	"${EQTL}" \
	"regions/${NAME}.chr${i}.regions.txt" \
	"Ashkenazi_LD_matrices/" \
	"${NAME}" \
	367

done

echo "We out!"
