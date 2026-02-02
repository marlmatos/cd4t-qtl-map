#!/bin/bash
#SBATCH --job-name=atac_workflow           # Job name
#SBATCH --nodes=1                          # node count
#SBATCH --partition=cpu                   # Partition Name
#SBATCH --mem=100G                          # Memory per node (increased for merge)
#SBATCH --cpus-per-task=8                  # CPU cores for the task
#SBATCH --time=2-00:00:00                    # Time limit (HH:MM:SS)
#SBATCH --output=logs/atac_prep_chrombnet_%j.log # Standard output log
#SBATCH --error=logs/atac_prep_chrombnet_%j.err  # Standard error log
#SBATCH --mail-type=END,FAIL              # Mail events
#SBATCH --mail-user=mmatos@nygenome.org   # Where to send mail

# Load necessary modules
module load singularity
module load samtools
module load bedtools

# Create logs directory if it doesn't exist
mkdir -p logs

# Create output directory
output_dir="/gchm/cd4_chrombpnet/data/inputs"
mkdir -p "$output_dir"

# Input directory where BAM files are located
bam_dir="/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/filtered_bams"

# Path to the text file containing sample names
samples_file="/gchm/cd4_chrombpnet/data/atac_high_read_depth_samples.txt"

# Reference files
blacklist_file="/gchm/resources/genome/hg38_gencode_raw/hg38-blacklist.v3.bed" 
chrom_sizes_file="/gchm/resources/genome/hg38.p14.chrom.sizes.fmtd" 

# Path to the singularity image
singularity_img="/gchm/.singularity/singularity_images_nextflow/depot.galaxyproject.org-singularity-macs3-3.0.1--py310h1af8fb7_3.img"

# Set genome size (human genome)
genome_size="2.86e9"

# Parameters for peak calling
p_value="0.01"
ext_size="183"

echo "Starting ATAC-seq prep workflow for ChrmBPNet"
echo "=================================="

# Step 1: Create a list of BAM files from the sample file
echo "Creating BAM file list from samples in $samples_file"
rm -f atac_high_read_depth_samples.txt

while read -r sample; do
    bam_file="${bam_dir}/${sample}.bfilt.orpham.sorted.bam"
    if [ -f "$bam_file" ]; then
        echo "$bam_file" >> atac_high_read_depth_samples.txt
    else
        echo "Warning: BAM file not found for sample $sample: $bam_file"
    fi
done < "$samples_file"

num_files=$(wc -l < atac_high_read_depth_samples.txt)
echo "Found $num_files BAM files for merging"

if [ "$num_files" -eq 0 ]; then
    echo "Error: No BAM files found. Exiting."
    exit 1
fi

# Step 2: Merge BAM files
echo "Merging BAM files..."
merged_unsorted="${output_dir}/merged_cd4_atac_unsorted.bam"
merged_unsorted_filtered="${output_dir}/merged_cd4_atac_unsorted_autosomes.bam"
merged_bam="${output_dir}/merged_cd4_atac.bam"

# Create a space-separated list of BAM files
bam_list=$(cat atac_high_read_depth_samples.txt | tr '\n' ' ')

echo "Running: samtools merge -f $merged_unsorted $bam_list"
samtools merge -f "$merged_unsorted" $bam_list

# Filter BAM file to keep only autosomal chromosomes (chr1-chr22)
echo "Filtering BAM to keep only autosomal chromosomes..."
samtools view -h "$merged_unsorted" | awk '{if($0 ~ /^@/ || $3 ~ /^chr[1-9]$/ || $3 ~ /^chr1[0-9]$/ || $3 ~ /^chr2[0-2]$/){print $0}}' | samtools view -Shb - > "$merged_unsorted_filtered"

# Sort the merged BAM file
echo "Sorting merged BAM file..."
samtools sort -@8 "$merged_unsorted_filtered" -o "$merged_bam"

# Index the merged BAM file
echo "Indexing merged BAM file..."
samtools index "$merged_bam"

# Remove the unsorted BAM file to save space
rm "$merged_unsorted"

echo "Merged BAM created: $merged_bam"

# Filter chromosome sizes to include only autosomes (chr1-chr22)
echo "Filtering chromosome sizes file to include only autosomes..."
filtered_chrom_sizes="${output_dir}/hg38.autosomes.chrom.sizes"
grep -E '^chr([1-9]|1[0-9]|2[0-2])\s' "$chrom_sizes_file" | sort -k1,1V > "$filtered_chrom_sizes"
# Step 3: Run MACS3 peak calling on the merged BAM
echo "Running MACS3 peak calling on merged BAM..."
peak_dir="${output_dir}/peaks"
mkdir -p "$peak_dir"

echo "Using MACS3 singularity container: $singularity_img"

# Run MACS3 peak calling using singularity
singularity exec \
    --bind "${output_dir}:${output_dir}" \
    --bind "${peak_dir}:${peak_dir}" \
    "$singularity_img" \
    macs3 callpeak \
    -t "$merged_bam" \
    -f BAMPE \
    -g "$genome_size" \
    -n "merged_cd4_samples" \
    --outdir "$peak_dir" \
    --pvalue "$p_value" \
    --nomodel \
    --extsize "$ext_size" \
    --call-summits \
    2> "${peak_dir}/merged_cd4_samples_macs3.log"

# Check if peak calling was successful
if [ ! -f "${peak_dir}/merged_cd4_samples_peaks.narrowPeak" ]; then
    echo "Error: MACS3 peak calling failed. No peaks file found."
    exit 1
fi

echo "MACS3 peak calling completed"

# Step 4: Filter and prepare blacklist
echo "Filtering blacklist to include only autosomal regions..."
filtered_blacklist="${peak_dir}/blacklist_autosomes.bed"
grep -E '^chr([1-9]|1[0-9]|2[0-2])\s' "$blacklist_file" > "$filtered_blacklist"

# Step 5: Filter peaks against blacklist
echo "Filtering peaks against blacklist regions..."

# Ensure filtered blacklist file exists
if [ ! -f "$filtered_blacklist" ]; then
    echo "Error: Filtered blacklist file not found: $filtered_blacklist"
    exit 1
fi

# Ensure filtered chromosome sizes file exists
if [ ! -f "$filtered_chrom_sizes" ]; then
    echo "Error: Filtered chromosome sizes file not found: $filtered_chrom_sizes"
    exit 1
fi

# Extend blacklist regions by 1057 bp on both sides
blacklist_extended="${peak_dir}/blacklist_extended.bed"
bedtools slop -i "$filtered_blacklist" -g "$filtered_chrom_sizes" -b 1057 > "$blacklist_extended"

# Remove peaks that intersect with extended blacklist
filtered_peaks="${peak_dir}/merged_cd4_samples_peaks_no_blacklist.bed"
bedtools intersect -v -a "${peak_dir}/merged_cd4_samples_peaks.narrowPeak" -b "$blacklist_extended" > "$filtered_peaks"

# Sort the filtered peaks by chromosome and position
sorted_filtered_peaks="${peak_dir}/merged_cd4_samples_peaks_no_blacklist.sorted.bed"
sort -k1,1V -k2,2n "$filtered_peaks" > "$sorted_filtered_peaks"

# Count peaks before and after filtering
total_peaks=$(wc -l < "${peak_dir}/merged_cd4_samples_peaks.narrowPeak")
filtered_total=$(wc -l < "$filtered_peaks")
sorted_total=$(wc -l < "$sorted_filtered_peaks")
removed=$((total_peaks - filtered_total))

echo "Total peaks before filtering: $total_peaks"
echo "Total peaks after filtering: $filtered_total"
echo "Total peaks after sorting: $sorted_total"
echo "Removed peaks: $removed"

# Clean up intermediate files
rm "$blacklist_extended"

echo "Analysis complete. Final peaks file: $sorted_filtered_peaks"