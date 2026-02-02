#!/bin/bash
#SBATCH --job-name=2.splitfolds              # Job name
#SBATCH --partition=gpu                      # Partition Name
#SBATCH --mem=8G
#SBATCH --gres=gpu
#SBATCH --output=logs/2.splitfolds_%A.log               # Standard output and error log
#SBATCH --error=logs/2.splitfolds_%A.err
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mmatos@nygenome.org      # Where to send mail 

# Load required modules
module load chrombpnet
module load bedtools

# Create output directories
mkdir -p /gchm/cd4_chrombpnet/data/splits/
mkdir -p logs

# Path to input files
peaks_file="/gchm/cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed"
genome_file="/gchm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa"
chrom_sizes="/gchm/cd4_chrombpnet/data/inputs/hg38.autosomes.chrom.sizes" 
blacklist_file="/gchm/cd4_chrombpnet/data/inputs/peaks/blacklist_autosomes.bed" 

# Predefined test and validation chromosomes for each fold
test_fold0=("chr1" "chr3" "chr6" "chr14")
test_fold1=("chr2" "chr8" "chr9" "chr16")
test_fold2=("chr4" "chr11" "chr12" "chr15")
test_fold3=("chr5" "chr10" "chr18" "chr20" "chr22")
test_fold4=("chr7" "chr13" "chr17" "chr19" "chr21")

val_fold0=("chr8" "chr20")
val_fold1=("chr12" "chr17")
val_fold2=("chr22" "chr7")
val_fold3=("chr6" "chr21")
val_fold4=("chr10" "chr18")

# Base directory for outputs
output_dir="/gchm/cd4_chrombpnet/data/splits"
output_dir2="/gchm/cd4_chrombpnet/data/folds"

# Loop through each fold (0-4)
for fold_num in {0..4}; do
    echo "Processing fold ${fold_num}..."
    
    # Get the test and validation chromosomes for this fold
    # Convert array to space-separated string for the command arguments
    test_var="test_fold${fold_num}[@]"
    val_var="val_fold${fold_num}[@]"
    
    test_chrs="${!test_var}"
    valid_chrs="${!val_var}"
    
    echo "Test chromosomes: $test_chrs"
    echo "Validation chromosomes: $valid_chrs"
    
    # Define the output path correctly
    fold_output="${output_dir}/fold_${fold_num}"
    
    # Step 1: Generate chromosome splits
    echo "Generating chromosome splits for fold ${fold_num}..."
    chrombpnet prep splits \
        -c "$chrom_sizes" \
        -tcr $test_chrs \
        -vcr $valid_chrs \
        -op "$fold_output"
    
    # Check if the splits generation was successful
    if [ ! -f "${fold_output}.json" ]; then
        echo "Error: Failed to generate splits for fold ${fold_num}"
        continue
    fi
    
    # Step 2: Prepare nonpeaks using the generated splits file
    echo "Preparing nonpeaks for fold ${fold_num}..."
    chrombpnet prep nonpeaks \
        -g "$genome_file" \
        -p "$peaks_file" \
        -c "$chrom_sizes" \
        -fl "${fold_output}.json" \
        -br "$blacklist_file" \
        -o "${output_dir2}/f${fold_num}_output"
    
    echo "Completed processing for fold ${fold_num}"
    echo "----------------------------------------"
done

echo "All folds processed successfully"
