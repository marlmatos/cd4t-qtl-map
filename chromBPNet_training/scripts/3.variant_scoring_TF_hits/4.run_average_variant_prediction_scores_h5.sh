#!/bin/bash
#SBATCH --job-name=averg_allele_prediction           # Job name
#SBATCH --nodes=1                          # node count
#SBATCH --partition=cpu                   # Partition Name
#SBATCH --mem=100G                          # Memory per node (increased for merge)
#SBATCH --cpus-per-task=8                  # CPU cores for the task
#SBATCH --time=02:00:00                    # Time limit (HH:MM:SS)


# Load necessary modules
source activate renviron_ne1

OUT_PATH="~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5"
SCRIPT_PY="~/cd4t-qtl-map/chromBPNet_training/scripts/3.variant_scoring_TF_hits/4.average_variant_prediction_scores_h5.py"


python "$SCRIPT_PY" \
  --fold-glob "~/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold_perchrom/fold_*" \
  --per-fold-file-glob "**/cd4_tcells_AJ_common_variants.chr*.variant_predictions.h5" \
  --tsv-path "~/cd4_chrombpnet/scripts/3.variant_scoring_TF_hits/input/converted_variants_inpeaks.tsv" \
  --out "$OUT_PATH" \
  --keep-keys "observed/allele1_pred_counts,observed/allele1_pred_profiles,observed/allele2_pred_counts,observed/allele2_pred_profiles"