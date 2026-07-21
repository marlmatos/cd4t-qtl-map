#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

# 01. Define regions containing GWS SNPs
# Provide eQTL sumstats path
# eQTL name

numbers=($(seq 1 22))

for i in "${numbers[@]}";do

echo "Running chr${i}"

module purge
module load R/4.4.1

## NAME - name of eQTL analysis
## INDEP - path to conditionally independent variants by chr (should be a csv file)
## EQTL - path to eQTL summary statistics by chromosome (parquet file)

TFS=(ATF1 BCL11B CREB1 CTCF FOXO1 IKZF3 IRF1 IRF7 IRF9 KLF4 KLF5 NFATC3 NFYA NRF1 RELA REST RFX3 RUNX2 SP2 STAT1 STAT2 TCF7 THAP11 ZNF76)
DIR=$(pwd)
NAME="ATF1"
INDEP="/gchm/cd4_QTL_analysis/04_Run_interaction_eQTLs/interaction_TF_disrupted_chromBPNet/results/interaction_chrombpnet_tfs/TF_activity_interaction_chrombpnet_tfs/ATF1/ATF1_chr${i}.cis_qtl_top_assoc.txt.gz"
EQTL="/gchm/cd4_QTL_analysis/04_Run_interaction_eQTLs/interaction_TF_disrupted_chromBPNet/results/interaction_chrombpnet_tfs/TF_activity_interaction_chrombpnet_tfs/ATF1/ATF1_chr${i}.cis_qtl_pairs.chr${i}.parquet"

Rscript 01.Define_finemap_regions.eQTL.R \
	"${INDEP}" \
	"${NAME}" \
	"${DIR}"
	

# Calculate LD for regions containing significant SNPs
sbatch --wait 02.Generate_LD.sh "regions/${NAME}.chr${i}.regions.txt" \
			 "/gchm/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/per_chr_plink_files/chr${i}_ashkenazi.367.AF1.QC.BA.king2.hwe.annot" \
			 "Ashkenazi"

# Wait for job to finish
wait

# Finemap regions
sbatch 03.submit_finemap.sh \
	"${EQTL}" \
	"regions/${NAME}.chr${i}.regions.txt" \
	"Ashkenazi_LD_matrices/" \
	"${NAME}" \
	367

done

echo "Done!"
