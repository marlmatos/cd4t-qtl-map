#!/bin/bash
#SBATCH --job-name=var_concat_summary  # Job name
#SBATCH --partition=gpu                # Partition Name
#SBATCH --mem=10G
#SBATCH --gres=gpu
#SBATCH --time=24:00:00                # total run time limit (HH:MM:SS)
#SBATCH --output=logs/var_summary_%A.log  # Standard output with array job ID
#SBATCH --error=logs/var_summary_%A.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL           # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org  # Where to send mail

# Load modules for singularity
module purge
source /usr/share/lmod/lmod/init/bash
module load singularity

# Define full home path to replace tilde
HOME_DIR=/gchm

# Define singularity image
singularity_image=$HOME_DIR/packages/chrombpnet_latest.sif

# Define directories
PERFOLD_PERCHROM_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold_perchrom
PERFOLD_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold
OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores
MODIFIED_HELPERS_DIR="/gchm/packages/variant-scorer/src/"

# Create output directories
mkdir -p $PERFOLD_DIR
mkdir -p $OUTPUT_DIR

# Step 1: Concatenate chromosome results for each fold
echo "Concatenating chromosome results for each fold..."

for FOLD in {0..4}; do
    echo "Processing fold_${FOLD}..."
    FOLD_DIR=$PERFOLD_PERCHROM_DIR/fold_${FOLD}
    FOLD_OUTPUT_DIR=$PERFOLD_DIR/fold_${FOLD}
    mkdir -p $FOLD_OUTPUT_DIR
    
    # Initialize an empty file for concatenation
    CONCAT_FILE=$FOLD_OUTPUT_DIR/cd4_tcells_AJ_common_variants_variant_scores.tsv
    
    # Get the header from the first chromosome file (chr1)
    CHR1_FILE=$(find $FOLD_DIR -name "cd4_tcells_AJ_common_variants.chr1.variant_scores.tsv" | head -1)
    
    if [ -f "$CHR1_FILE" ]; then
        # Extract header
        head -1 $CHR1_FILE > $CONCAT_FILE
        
        # Concatenate all chromosome files (skipping headers)
        for CHR_FILE in $FOLD_DIR/cd4_tcells_AJ_common_variants.chr*.variant_scores.tsv; do
            if [ -f "$CHR_FILE" ]; then
                echo "  Adding $(basename $CHR_FILE)..."
                tail -n +2 $CHR_FILE >> $CONCAT_FILE
            fi
        done
        
        echo "  Created $CONCAT_FILE"
    else
        echo "  Warning: Could not find chr1 file for fold_${FOLD}"
    fi
done

# Step 2: Run the summary across folds
echo "Running summary across folds..."

# Define script paths
summary_script=$HOME_DIR/packages/variant-scorer/src/variant_summary_across_folds.py
annotation_script=$HOME_DIR/packages/variant-scorer/src/variant_annotation.py

# Define other paths
peaks=$HOME_DIR/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed
genes=$HOME_DIR/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed
out_prefix=$OUTPUT_DIR/average_cd4_tcells_AJ_common_variants

# Create a list of concatenated score files for all folds
SCORE_FILES=""
for FOLD in {0..4}; do
    SCORE_FILE=$PERFOLD_DIR/fold_${FOLD}/cd4_tcells_AJ_common_variants_variant_scores.tsv
    if [ -f "$SCORE_FILE" ]; then
        SCORE_FILES="$SCORE_FILES $SCORE_FILE"
    fi
done

echo "Using score files: $SCORE_FILES"

# Run the summary script using singularity
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --env MPLCONFIGDIR=/tmp \
  --env PYTHONPATH="$MODIFIED_HELPERS_DIR:$PYTHONPATH" \
  --no-home \
  -B $HOME_DIR/packages \
  -B $HOME_DIR/cd4_chrombpnet \
  -B $HOME_DIR/resources \
  -B /gchm/ \
  -B $OUTPUT_DIR \
  -B $PERFOLD_DIR \
  $singularity_image python $summary_script \
  --score_dir $PERFOLD_DIR \
  --score_list $SCORE_FILES \
  --out_prefix $out_prefix \
  --schema 'chrombpnet'

echo "Done!"