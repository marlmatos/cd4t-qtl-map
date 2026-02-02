#!/bin/bash
#SBATCH --job-name=hitcaller       # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=8G
#SBATCH --gres=gpu
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output=logs/finemo.log_%A_%a.log  # Standard output with array job ID
#SBATCH --error=logs/finemo.log_%A_%a.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail

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

# Continue with your original script
HOME_DIR=/gpfs/commons/home/mmatos
OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs/motif_instances
mkdir -p $OUTPUT_DIR

shap_data=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged/averaged_folds_cd4_tcells.counts_scores.h5
modisco=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs/from_count_contributions/tfmodisco_motifs_count_contributions.h5

narrowpeaks="/gpfs/commons/home/mmatos/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed"
hits_file="$OUTPUT_DIR/hits.tsv"

echo "Extracting Regions"
finemo extract-regions-chrombpnet-h5 -c $shap_data -o ${OUTPUT_DIR}/regions_input -w 400

echo "Calling Instances"
finemo call-hits -r  ${OUTPUT_DIR}/regions_input.npz -m $modisco -o ${OUTPUT_DIR} -p $narrowpeaks

echo "Constructing Report"
finemo report --regions ${OUTPUT_DIR}/regions_input.npz --hits $hits_file --out-dir ${OUTPUT_DIR} --modisco-h5 $modisco --modisco-region-width 400 -p $narrowpeaks

echo "Finished"