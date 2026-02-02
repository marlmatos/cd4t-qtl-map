#!/bin/bash
#SBATCH --job-name=prediction_scores         # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=100G
#SBATCH --gres=gpu
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output=logs/prediction_scores_%A_%a.log  # Standard output with array job ID
#SBATCH --error=logs/prediction_scores_%A_%a.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail
#SBATCH --array=0-4                         # Array range - adjust based on number of folds/splits

#this step is pretty quick 

source /usr/share/lmod/lmod/init/bash

#  full home path 
HOME_DIR=/gchm

#  the singularity image path
singularity_image=$HOME_DIR/packages/chrombpnet_latest.sif  # Adjust this path to where your image is stored

# Get the current array task ID
FOLD_ID=${SLURM_ARRAY_TASK_ID}

#  output directory first so we can use it for other path definitions
OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/prediction_scores_bw/fold_${FOLD_ID}/

genome=$HOME_DIR/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa
model=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet.h5
model_nobias=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet_nobias.h5
chromsizes=$HOME_DIR/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes
peaks=$HOME_DIR/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed

out_prefix=${OUTPUT_DIR}/cd4_tcells

# Create output directory
mkdir -p $OUTPUT_DIR

# Echo the parameters for debugging
echo "Processing fold: $FOLD_ID"
echo "Using ChromBPNet model: $model"
echo "Using bias-corrected model: $model_nobias"
echo "Output directory: $OUTPUT_DIR"
echo "Using singularity image: $singularity_image"

module load singularity

# Run the contribution scores command
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --no-home \
  -B $HOME_DIR/packages \
  -B $HOME_DIR/cd4_chrombpnet \
  -B $HOME_DIR/resources \
  $singularity_image \
  chrombpnet pred_bw \
  -g $genome \
  -cmb $model_nobias \
  -r $peaks \
  -c $chromsizes \
  -op $out_prefix 