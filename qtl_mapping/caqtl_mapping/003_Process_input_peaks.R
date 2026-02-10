#------------------------
# Description: This script prepare input peak file for tensorQTL
#
# Author: Marliette Matos
# Date: 09/25/2024
# 
#--------------------------
#libraries
.libPaths("/gpfs/commons/home/mmatos/backup_R/x86_64-pc-linux-gnu-library")
library(data.table)
library(dplyr)
library(tidyr)
set.seed(2024)

# Read command-line arguments from SLURM job
args <- commandArgs(trailingOnly = TRUE)
filtered_table <- args[1]  # This is the table file passed from the SLURM array
filtered_table="filtered_qsmooth_norm_cpm.txt"
# Directories
root_dir <- "/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/"
annot_file <- "/gpfs/commons/groups/lappalainen_lab/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/merged_library/peak_calling/MACS3/BAMPE/peaks_102024/cd4_atac_padded_summits_peaks.bed" #selecting this becuase already has the correct format

# Make output directory for the specific table
out_dir <- paste0(root_dir, "003_inputs/", gsub(".txt", "", filtered_table), "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# Get Peak quantification file (VST or CPM + TMM)
peaks_path <- paste0(root_dir, "001_peaks/", filtered_table)
print(paste("Processing peaks from:", peaks_path))

peaks <- fread(peaks_path, header = TRUE)
peaks <- as.data.frame(peaks)
colnames(peaks) <- paste0("g", colnames(peaks))
rownames(peaks) <- peaks[[1]] # set the first column as the rownames
peaks <- peaks[ , -1] # remove the first column
peaks=as.data.frame(t(peaks))

peaks[1:5,1:5]

print("Remove individuals with NA")
peaks <- mutate_all(peaks, function(x) as.numeric(as.character(x)))
na_idx <- which(is.na(rowMeans((peaks))))
if(length(na_idx) > 0){
  peaks <- peaks[-na_idx, ]
}
                    
dim(peaks)
print("Remove peaks with all zeros")
na_idx <- which((colSums((peaks)) == 0))
if(length(na_idx) > 0){
  peaks <- peaks[ , -na_idx]
}
dim(peaks)

#transform the counts back so the the peaks are each row
peaks <- t(peaks)

# Remove peaks with pi0 > 0.9
pi0 = rowSums(peaks==0)/ncol(peaks)
peaks = peaks[which(pi0 <= 0.9), ]
dim(peaks)

                    
# log(x+1) the cpm counts and standardization                    
peaks <- log2(peaks + 1)


# Convert matrix to BED format

# Get the peak coordinates
# Read the first few columns from the annotation file
annot = fread(annot_file)
colnames(annot)= c('chr', 'start','end', 'peak_name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak', 'peak_summit')
                    
peak_coords1 = annot[annot$peak_name %in% rownames(peaks), c(1, 2, 3, 4)]
write.table(peak_coords1, paste0(out_dir, "cd4_atac_processed_peaks_coordinates.bed"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# Use the peak names to get the coordinates of the peaks from the annotation file
peak_coords = annot[annot$peak_name %in% rownames(peaks), c(4, 1, 11)]

# Format the peaks file to have columns chr, start, end, gene_name where the end is just start + 1
peak_coords$start = as.numeric(peak_coords$peak_summit)
peak_coords$end = as.numeric(peak_coords$peak_summit) + 1

# Reorder the columns
peak_coords = peak_coords[, c(1, 2, 4, 5)]

colnames(peak_coords)= c('peak_name', 'chr','start', 'end')  # Add gene_name column from rownames

data = cbind(peak_coords[, c("chr", "start", "end", "peak_name")], peaks[,-1])
data = data[order(data$chr, data$start), ]
colnames(data)[1] = "#chr"


# Remove NAs
data = na.omit(data)
write.table(data, paste0(out_dir, gsub(".txt", "", filtered_table), "_cd4_atac_processed_peaks.bed"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
cat("Done!", fill = TRUE)

sessionInfo()