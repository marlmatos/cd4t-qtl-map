#!/bin/bash
#SBATCH --job-name=chrombpnet_array         # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=100G
#SBATCH --gres=gpu
#SBATCH --time=24:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --array=0-4                         # Array range - adjust based on number of folds/splits

source /usr/share/lmod/lmod/init/bash
module load singularity


# singularity image
singularity_image=$HOME/packages/chrombpnet_latest.sif

# Get the current array task ID
FOLD_ID=${SLURM_ARRAY_TASK_ID}
# output directory with array index
OUTPUT_DIR="$HOME/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/"

# Create output directory
mkdir -p $OUTPUT_DIR

# fold-specific paths
FOLD_FILE="$HOME/cd4_chrombpnet/data/splits/fold_${FOLD_ID}.json"
NEGATIVE_FILE="$HOME/cd4_chrombpnet/data/folds/f${FOLD_ID}_output_negatives.bed"
# Bias model stays constant
BIAS_MODEL="$HOME/cd4_chrombpnet/bias_model_f0_b7/models/cd4_bias.h5"

# Echo the parameters for debugging
echo "Processing fold: $FOLD_ID"
echo "Using fold file: $FOLD_FILE"
echo "Using negative file: $NEGATIVE_FILE"
echo "Using bias model: $BIAS_MODEL"
echo "Output directory: $OUTPUT_DIR"

# Run singularity with environment isolation
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --no-home \
  -B $HOME/packages \
  -B $HOME/cd4_chrombpnet \
  -B $HOME/resources \
  $singularity_image \
  chrombpnet pipeline \
        -ibam $HOME/cd4_chrombpnet/data/inputs/merged_cd4_atac.bam \
        -d "ATAC" \
        -g $HOME/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa \
        -c $HOME/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes \
        -p $HOME/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed \
        -n $NEGATIVE_FILE \
        -fl $FOLD_FILE \
        -b $BIAS_MODEL \
        -o $OUTPUT_DIR