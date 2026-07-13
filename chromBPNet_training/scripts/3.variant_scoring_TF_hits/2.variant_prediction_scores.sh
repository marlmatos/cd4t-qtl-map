#!/bin/bash
#SBATCH --job-name=variant_score         # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=100G
#SBATCH --gres=gpu
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --array=1-110%10  # 5 folds × 22 chromosomes = 110 total tasks, limit to 10 concurrent

source /usr/share/lmod/lmod/init/bash
module purge
module load singularity 

# Define chromosomes array
CHROMOSOMES=(
  "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10"
  "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" 
  "chr20" "chr21" "chr22"
)

# Number of chromosomes
NUM_CHROMS=22

# Calculate fold and chromosome from array task ID
# Integer division to get fold (0-based)
FOLD_IDX=$(( (SLURM_ARRAY_TASK_ID-1) / NUM_CHROMS ))
# Modulo to get chromosome index (0-based)
CHROM_IDX=$(( (SLURM_ARRAY_TASK_ID-1) % NUM_CHROMS ))

# Get fold ID (starting from 0) and chromosome name
FOLD_ID=$FOLD_IDX
CHROM=${CHROMOSOMES[$CHROM_IDX]}

# Echo for debugging
echo "Processing fold: $FOLD_ID, chromosome: $CHROM"

# Define paths based on fold
OUTPUT_DIR=$HOME/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/perfold_perchrom/fold_${FOLD_ID}/
mkdir -p $OUTPUT_DIR

list=$HOME/cd4_chrombpnet/scripts/3.variant_scoring_TF_hits/input/converted_variants_inpeaks.tsv
genome=$HOME/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa
model=$HOME/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet.h5
model_nobias=$HOME/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet_nobias.h5
chromsizes=$HOME/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes
peaks=$HOME/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed
out_prefix=${OUTPUT_DIR}/cd4_tcells_AJ_common_variants

MODIFIED_HELPERS_DIR="/gchm/packages/variant-scorer/src/"
variant_scoring_py=/gchm/packages/variant-scorer/src/variant_scoring.per_chrom.py
singularity_image=$HOME/packages/chrombpnet_latest.sif

# Run variant scoring for specific fold and chromosome
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --env MPLCONFIGDIR=/tmp \
  --env PYTHONPATH="$MODIFIED_HELPERS_DIR:$PYTHONPATH" \
  --no-home \
  -B $HOME/packages \
  -B $HOME/cd4_chrombpnet \
  -B $HOME/resources \
  -B /gchm/ \
  -B $OUTPUT_DIR \
  $singularity_image python $variant_scoring_py \
  -l $list \
  -g $genome \
  -m $model_nobias \
  -o $out_prefix \
  -s $chromsizes \
  -p $peaks \
  -sc 'chrombpnet' \
  -c $CHROM