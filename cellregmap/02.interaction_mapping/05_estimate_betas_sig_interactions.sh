#!/bin/bash

#SBATCH --job-name=beta_est                  # Job name
#SBATCH --nodes=1                # node count               
#SBATCH --partition=cpu               # Partition Name
#SBATCH --mem=200G
#SBATCH --cpus-per-task=4  # Request 4 CPU cores for the task
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output=logs/beta_est%A_%a.log               # Standard output and error log
#SBATCH --array=0-248
#SBATCH --error=logs/beta_est%A_%a.err

source activate cellregmap

## SLURM ARRAY ##
echo "$SLURM_ARRAY_TASK_ID"

res=34
fdr=5

python 05_estimate_betas_sig_interactions.py \
  --index ${SLURM_ARRAY_TASK_ID} \
  --res ${res} \
  --fdr ${fdr} \
  --input_files_dir "~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/" \
  --results_dir "~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/" \
  --metadata_file "~/metadata_harmonizing/results/post_WGS_QC/version_080923/revision_103023/cd4_scRNAseq_meta.v3.csv" \
  --mapping_file "~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_sample_mapping_res${res}.csv" \
  --contexts_file "~/cd4_CellRegMap/001_preprocessing/results_01312025/mofa_factors/mofa_trained/cd4_aging_filt_sce_mofa_expectations_model_factors_res${res}.csv" \
  --covariates_file "~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/Data/pseudobulk_metacells/pseusobulk_covariates_res${res}.csv" \
  --int_results_dir "~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/many_lead_snps_per_gene_res${res}_5factors/fdr_corrected" \
  --out_dir "~/cd4_CellRegMap/002_interaction_analysis/results/results_042426/many_lead_snps_per_gene_res${res}_5factors/beta_estimation/fdr${fdr}"


