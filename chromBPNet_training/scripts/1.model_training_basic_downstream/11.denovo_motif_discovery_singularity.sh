#!/bin/bash
#SBATCH --job-name=tfmodisco         # Job name
#SBATCH --partition=cpu                     # Partition Name
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00                     # total run time limit (HH:MM:SS)
#SBATCH -c 10
#SBATCH --output=logs/tfmodisco_%A.log  # Standard output with array job ID
#SBATCH --error=logs/tfmodisco_%A.err   # Error log with array job ID
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org     # Where to send mail

#please note that tfmodisco within chrombpnet will be deprecated. The reason I used the version within the singularity image is that I was having some dependency issue with numpy

source /usr/share/lmod/lmod/init/bash

# full home path
HOME_DIR=/gchm

# the singularity image path
singularity_image=$HOME_DIR/packages/chrombpnet_latest.sif  # Adjust this path to where your image is stored
# Set Numba to use 4 threads
export NUMBA_NUM_THREADS=8


# output directory first so we can use it for other path definitions
OUTPUT_DIR=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs/from_count_contributions
OUTPUT=$OUTPUT_DIR//tfmodisco_motifs_count_contributions.h5 
contribution_scores=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/contribution_scores_bw/averaged/averaged_folds_cd4_tcells.counts_scores.h5
model_nobias=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/fold_${FOLD_ID}/models/chrombpnet_nobias.h5
meme=$HOME_DIR/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs/motifs_database/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt 
# Create output directory
mkdir -p $OUTPUT_DIR/report/

# Echo the parameters for debugging
echo "Using averaged contribution scores from ChomBPNet model: $contribution_scores derived from $model_nobias"
echo "Output directory: $OUTPUT_DIR"
echo "Using singularity image: $singularity_image"

module load singularity

cd $OUTPUT_DIR

# Run the contribution scores command
#note that w=400 is the default
#to get a representative summary of all the TF-motifs learnt by the model, we recommend running tfmodisco-lite on contribution scores from all the peaks and with the -n argument set to 1000000 
# singularity exec --nv \
#   --env PYTHONNOUSERSITE=1 \
#   --env MPLCONFIGDIR=/tmp \
#   --no-home \
#   -B $HOME_DIR/packages \
#   -B $HOME_DIR/cd4_chrombpnet \
#   $singularity_image modisco motifs \
#   -i "$contribution_scores" \
#   -n 1000000 \
#   -w 400 \
#   --verbose \
#   -o $OUTPUT

#generate the report
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --env MPLCONFIGDIR=/tmp \
  --no-home \
  -B $HOME_DIR/packages \
  -B $HOME_DIR/cd4_chrombpnet \
  $singularity_image modisco report -i $OUTPUT -o report/ -s report/ -m $meme