#!/bin/bash
#SBATCH --job-name=3.bias_model                # Job name
#SBATCH --partition=gpu                      # Partition Name
#SBATCH --mem=100G
#SBATCH --gres=gpu
#SBATCH --time=2-00:00:00          # total run time limit (HH:MM:SS)


source /usr/share/lmod/lmod/init/bash
module load singularity

# Define singularity image
singularity_image=$HOME/packages/chrombpnet_latest.sif

mkdir -p $HOME/cd4_chrombpnet/bias_model_f0_b7

# Run singularity with environment isolation
singularity exec --nv \
  --env PYTHONNOUSERSITE=1 \
  --no-home \
  -B $HOME/packages \
  -B $HOME/cd4_chrombpnet \
  -B $HOME/resources \
  $singularity_image \
  chrombpnet bias pipeline \
  -ibam $HOME/cd4_chrombpnet/data/inputs/merged_cd4_atac.bam \
  -d "ATAC" \
  -g $HOME/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa \
  -c $HOME/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes \
  -p $HOME/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed \
  -n $HOME/cd4_chrombpnet/data/folds/f0_output_negatives.bed \
  -fl $HOME/cd4_chrombpnet/data/splits/fold_0.json \
  -b 0.7 \
  -o $HOME/cd4_chrombpnet/bias_model_f0_b7/ \
  -fp cd4
