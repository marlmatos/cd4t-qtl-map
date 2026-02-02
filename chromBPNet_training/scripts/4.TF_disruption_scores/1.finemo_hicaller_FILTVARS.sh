#!/bin/bash
#SBATCH --job-name=hitcaller       # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=50G
#SBATCH --gres=gpu
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output=logs/hitcaller.log_%A_%a.log  # Standard output with array job ID
#SBATCH --error=logs/hitcaller.log_%A_%a.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail
# Array range - adjust based on number of folds/splits

###
# This script is for identifyings instances of modisco motifs
###

# Clean the environment first
module purge
source /usr/share/lmod/lmod/init/bash

# Activate your finemo environment but with more isolation
source activate finemo

# Ensure Python only sees packages from the conda environment by setting PYTHONPATH
export PYTHONPATH=$CONDA_PREFIX/lib/python3.10/site-packages:$PYTHONPATH

# Set PATH to prioritize conda binaries
export PATH=$CONDA_PREFIX/bin:$PATH

HOME_DIR=/gpfs/commons/home/mmatos
OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/motif_hit_calls/motifs_hits_FILTVARS

mkdir -p $OUTPUT_DIR

# Original variant file
variants=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/cd4_top_cpbnet_variants_noheader.tsv

# Create a file with NO header - it's cleaner this way since the script will add column names
#variant_no_header=${OUTPUT_DIR}/variants_no_header.tsv
#tail -n +2 $variants > $variant_no_header

# Check the file
echo "Checking file format (first few lines):"
head -n 3 $variants

hitcaller_script=$HOME_DIR/packages/variant-scorer/src/hitcaller_variant.py
shap_data=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5
modisco=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs_hocomoco_jaspar_cisbp/model/tfmodisco_motifs_count_contributions.h5

python $hitcaller_script \
  --shap_data $shap_data \
  --input_type "h5" \
  --modisco_h5 $modisco \
  --variant_file $variants \
  --hits_per_loc 4 \
  --output_dir ${OUTPUT_DIR} \
  --alpha 0.6
  
  
  
