#------------------------
# Description: This script calculates PCs for the ATACseq peak quantification matrix for caQTL analysis
#
# Author: Marliette Matos
# Date: 10/14/2024
# 
#--------------------------
# Libraries
.libPaths("/gchm/backup_R/x86_64-pc-linux-gnu-library")
library(PCAForQTL)
library(data.table)
library(RNOmni)
library(dplyr)
set.seed(2024)

# Read command-line arguments from SLURM job
args <- commandArgs(trailingOnly = TRUE)
filtered_table <- args[1]  # This is the table file passed from the SLURM array

# Directories
root_dir <- "/gchm/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/"
metadata_path <- "/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/metadata/cd4_atac_metadata_scrna_union.csv"

filtered_table="filtered_qsmooth_norm_cpm.txt"
# Make output directory for the specific table
out_dir <- paste0(root_dir, "002_pca/", gsub(".txt", "", filtered_table), "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Get Peak quantification file (VST or CPM + TMM)
peaks_path <- paste0(root_dir, "001_peaks/", filtered_table)
print(paste("Processing peaks from:", peaks_path))

peaks <- fread(peaks_path, header = TRUE, )
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


# Read Metadata
metadata <- read.csv(metadata_path, row.names = 1)
metadata$WGS_sampleID <- paste0("g", metadata$WGS_sampleID)

# Filter metadata for common samples
metadata <- metadata[metadata$WGS_sampleID %in% rownames(peaks), ]
metadata[1:5, 1:5]

# Read Ashkenazi Ancestry Information
ashkenazi_ancestry <- read.csv("~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/scripts/AJ_ancestry.csv")
colnames(ashkenazi_ancestry) <- c("WGS_sampleID", "Parent_Directory", "AJ_ancestry_%")
ashkenazi_ancestry$WGS_sampleID <- paste0("g", ashkenazi_ancestry$WGS_sampleID)

#read other atac-seq qc metrics 
seq_metrics<- read.csv("/gcglm/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/atac_preprocessing_metrics_qc_combined.csv" )           

# Merge with Metadata
metadata <- merge(metadata, seq_metrics, by.x = "ATAC_Sample_Name", by.y="Library")                 
metadata <- left_join(metadata, ashkenazi_ancestry, by = "WGS_sampleID")

# Columns to keep in metadata
metadata <- metadata[, c("WGS_sampleID", "age", "Sex", "ATAC_Lane", "Fasting", "AJ_ancestry_%", "FRiP_Score", "Non_CD4_T_Proportion", "ATAC_Batch", "p_CD4_Naive", "p_CD4_Central_Memory", "p_CD4_Effector_Memory", "p_CD4_Cytotoxic", "p_CD4_Treg", "Weighted_Mean_GC_Content", "Mean_Base_Quality", "Mean_N_Content", "Total_Aligned_Filtered_Reads", "Picard_PCT_PF_READS", "Picard_PCT_PF_READS_ALIGNED", "Picard_PF_ALIGNED_BASES", "Picard_PF_HQ_ALIGNED_Q20_BASES", "Picard_PF_MISMATCH_RATE", "Picard_PF_HQ_ERROR_RATE", "Picard_MEAN_READ_LENGTH", "Picard_MEDIAN_READ_LENGTH", "Picard_STRAND_BALANCE", "Picard_PCT_READS_ALIGNED_IN_PAIRS")]

# Match metadata order with peaks
metadata <- metadata[match(rownames(peaks), metadata$WGS_sampleID), ]

# Remove duplicates
metadata <- metadata[!duplicated(metadata$WGS_sampleID), ]

# Set row names to WGS_sampleID
rownames(metadata) <- metadata$WGS_sampleID

write.csv(metadata, paste0(out_dir, "cd4_atac_metadata_filtered.csv"), row.names=FALSE)
                    
print("Remove peaks with pi0 > 0.9")
pi0 = rowSums(peaks==0)/ncol(peaks)
peaks = peaks[which(pi0 <= 0.9), ]
dim(peaks)

                    
# Log transform the CPM counts
print("Applying log2(x + 1) transformation to CPM peaks")

    # Apply log2(x + 1) to each column (since peaks are columns)
    peaks <- apply(peaks, 2, function(x) log2(x + 1))

                    
## PCA Analysis ##
print("Start the PCA estimation...")
tick1 <- Sys.time()

# Check for NA in peaks
print("Checking for NAs")
table(is.na(peaks))
na_idx <- which(is.na(colMeans((peaks))))
if(length(na_idx) > 0){
  peaks <- peaks[ , -na_idx]
}

                                        
# Run PCA
prcompResult <- prcomp(peaks, center = TRUE, scale = TRUE)
PCs <- prcompResult$x

tick2 <- Sys.time()
print(paste0("Running time for PCA estimation is ", as.numeric(tick2 - tick1, units = "mins"), " mins"))

## Select PCs based on different methods ##
# Elbow method
print("Selecting PCs with elbow method")
resultRunElbow <- PCAForQTL::runElbow(prcompResult = prcompResult)
K_elbow <- resultRunElbow

# BE method
print("Selecting PCs with BE method")
resultRunBE <- PCAForQTL::runBE(peaks, B = 20, alpha = 0.05)
K_BE <- resultRunBE$numOfPCsChosen

# Plot and save scree plot
PCAForQTL::makeScreePlot(prcompResult, labels = c("Elbow", "BE"), values = c(K_elbow, K_BE), titleText = "CD4 T Cells Narrow Peaks")
ggplot2::ggsave(paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_peaks_PCs.pdf'), width = 16, height = 11, units = "cm")

# Save Elbow method results
factors_df_elbow <- cbind(sampleid = rownames(PCs), metadata, PCs[, 1:K_elbow])
write.table(factors_df_elbow, paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_peaks_elbow_pca_', K_elbow, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)
write.csv(factors_df_elbow, file = paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_peaks_elbow_pca_', K_elbow, '.csv'), row.names = FALSE)

# Save BE method results
factors_df_BE <- cbind(sampleid = rownames(PCs), metadata, PCs[, 1:K_BE])
write.table(factors_df_BE, paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_peaks_BE_pca_', K_BE, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE)
write.csv(factors_df_BE, file = paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_peaks_BE_pca_', K_BE, '.csv'), row.names = FALSE)

end <- Sys.time()
print(paste0("Running time for whole script is ", as.numeric(end - tick1, units = "mins"), " mins"))

sessionInfo()