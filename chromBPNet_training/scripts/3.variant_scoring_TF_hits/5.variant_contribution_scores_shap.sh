#!/bin/bash
#SBATCH --job-name=variant_shap       # Job name
#SBATCH --partition=gpu                     # Partition Name
#SBATCH --mem=50G
#SBATCH --gres=gpu
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH --output=logs/variant_shap.log_%A_%a.log  # Standard output with array job ID
#SBATCH --error=logs/variant_shap.log_%A_%a.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail
#SBATCH --array=0-4                         # Array range - adjust based on number of folds/splits

module purge
source /usr/share/lmod/lmod/init/bash
module load singularity 

# Ensure Python only sees packages from the conda environment by setting PYTHONPATH
export PYTHONPATH=$CONDA_PREFIX/lib/python3.10/site-packages:$PYTHONPATH

# Set PATH to prioritize conda binaries
export PATH=$CONDA_PREFIX/bin:$PATH


HOME_DIR=/gchm
singularity_image=$HOME_DIR/packages/chrombpnet_latest.sif
FOLD_ID=${SLURM_ARRAY_TASK_ID}

OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/perfold/fold_${FOLD_ID}/
mkdir -p $OUTPUT_DIR

# Original variant file
variants=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/cd4_top_cpbnet_variants_noheader.tsv
genome=$HOME_DIR/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa
model_nobias=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet_nobias.h5
chromsizes=$HOME_DIR/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes
out_prefix=cd4_tcells_AJ_common_variants
shap_script=$HOME_DIR/packages/variant-scorer/src/variant_shap.py
# Define the directory containing your modified helpers.py
HELPERS_DIR="$HOME_DIR/packages/variant-scorer/src/utils"

echo "Processing fold: $FOLD_ID"

# First check if the modified helpers file exists
if [ ! -f "$HELPERS_DIR/helpers.py" ]; then
  echo "Error: helpers file not found at $HELPERS_DIR/helpers.py"
  exit 1
fi


# Now run the singularity container with PYTHONPATH set to include your modified helpers directory
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --env MPLCONFIGDIR=/tmp \
  --env PYTHONPATH="$HELPERS_DIR:$PYTHONPATH" \
  --no-home \
  -B $HOME_DIR/packages \
  -B $HOME_DIR/cd4_chrombpnet \
  -B $HOME_DIR/resources \
  -B /gchm/ \
  -B $OUTPUT_DIR \
  $singularity_image python $shap_script \
  -l $variants \
  -g $genome \
  -s $chromsizes \
  -m $model_nobias \
  -o ${OUTPUT_DIR}/${out_prefix}.shap \
  -sc chrombpnet