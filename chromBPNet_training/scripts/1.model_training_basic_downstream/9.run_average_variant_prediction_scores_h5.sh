#!/bin/bash
#SBATCH --job-name=averg_allele_prediction           # Job name
#SBATCH --nodes=1                          # node count
#SBATCH --partition=cpu                   # Partition Name
#SBATCH --mem=100G                          # Memory per node (increased for merge)
#SBATCH --cpus-per-task=8                  # CPU cores for the task
#SBATCH --time=02:00:00                    # Time limit (HH:MM:SS)
#SBATCH --output=logs/averg_allele_prediction_h5_%j.log # Standard output log
#SBATCH --error=logs/averg_allele_prediction_h5_%j.err  # Standard error log
#SBATCH --mail-type=END,FAIL              # Mail events
#SBATCH --mail-user=mmatos@nygenome.org   # Where to send mail

# Load necessary modules
source activate renviron_ne1

OUT_PATH="/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_variant_prediction_scores.h5"
SCRIPT_PY="/gchm/cd4_chrombpnet/scripts/1.model_training_basic_downstream/9.average_variant_prediction_scores_h5.py"

# Optional: print what will be processed (no write)
python "$SCRIPT_PY" \
  --fold-glob "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold_perchrom/fold_*" \
  --per-fold-file-glob "**/cd4_tcells_AJ_common_variants.chr*.variant_predictions.h5" \
  --out "$OUT_PATH" \
  --dataset-prefix "observed" \
  --keep-keys "observed/allele1_pred_counts,observed/allele1_pred_profiles,observed/allele2_pred_counts,observed/allele2_pred_profiles" \
  --dry-run

# Real run (remove --dry-run)
python "$SCRIPT_PY" \
  --fold-glob "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold_perchrom/fold_*" \
  --per-fold-file-glob "**/cd4_tcells_AJ_common_variants.chr*.variant_predictions.h5" \
  --out "$OUT_PATH" \
  --dataset-prefix "observed" \
  --keep-keys "observed/allele1_pred_counts,observed/allele1_pred_profiles,observed/allele2_pred_counts,observed/allele2_pred_profiles"

