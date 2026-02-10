#-------------------------------------------------------
# Create covariates file for tensorQTL
# Author: Winona Oliveros
# Adpated by: Marliette Matos
#------------------------------------------------------

#Set up environment-------------------------------------------------------
#libraries
library(dplyr)
library(tidyr)
library(data.table)
library(Hmisc)
library(corrplot)
library(ggplot2)
library(RColorBrewer)
set.seed(2024)

# Read command-line arguments from SLURM job
args <- commandArgs(trailingOnly = TRUE)
filtered_table <- args[1]  # This is the table file passed from the SLURM array
filtered_table="filtered_qsmooth_norm_cpm.txt"

#directories
root_dir="/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/"

# Make output directory for the specific table
out_dir <- paste0(root_dir, "005_covariates/", gsub(".txt", "", filtered_table), "/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#Prepare covariates file for tensorQTL-----------------------

## process Genotypes PCs
genotype <- read.table(paste0(root_dir,'004_genotypes/plink/CD4_all_chr_ashkenazi.362.AF1.QC.BA.king2.hwe.ld.eigenvec'), header = F) ### genotype PC are the same for all celltypes

genotype <- genotype %>% mutate(sampleid = paste0("g", V2 ))
colnames(genotype)[colnames(genotype) == "V2"] <- "sampleid"
genotype <- genotype[, !colnames(genotype) %in% "V1"]
#-----------subset to adequate number of genotype PCS
genotype <- genotype[,c("sampleid","V3","V4","V5","V6","V7","V8", "V9", "V10","V11","V12","V13", "V14", "V15", "V16", "V17")]
genotype$sampleid <- paste0("g", genotype$sampleid)


# Replace 'V' with 'Genotype_PC' in each column name that starts with 'V'
# Modify column names that start with 'V' and increment their numbers
new_names <- sapply(colnames(genotype), function(x) {
  if(grepl("^V\\d+", x)) {
    num <- as.numeric(sub("^V", "", x))  # Extract the number part after 'V'
    return(sprintf("Genotype_PC%d", num - 2))  # Adjust the number
  } else {
    return(x)  # Return the name unchanged if it does not start with 'V'
  }
})

# Assign the new names back to the dataframe
colnames(genotype) <- new_names


#----------- Read in ATAC-seq metadata for Experimental Covariates 
peak_metadata="cd4_atac_metadata_filtered.csv"
metadata<-read.csv2(paste0(root_dir, "002_pca/",  gsub(".txt", "", filtered_table),"/" , peak_metadata), sep=",")

# process chromatin accessiility feature PCs
# Construct the base path without the number of PCs
base_path <- paste0(root_dir, "002_pca/", gsub(".txt", "", filtered_table),"/", gsub(".txt", "", filtered_table), "_atac_peaks_elbow_pca_")

# List files that match the pattern (find all files that start with base_path and end with .csv)
pcs_files <- list.files(dirname(base_path), pattern = paste0(basename(base_path), ".*\\.csv"), full.names = TRUE)

# Check if any files were found
if (length(pcs_files) > 0) {
    # If there's more than one file, use the first one (or decide based on other criteria)
    pcs_path <- pcs_files[1]  # Choose the first file, or add logic to select the best match

    # Read the CSV file
    PCS_1 <- read.csv(pcs_path, sep = ",")
    print(paste("Successfully loaded:", pcs_path))
} else {
    stop("No matching file found.")
}

#covariates_PCs <- PCS_1[,c("sampleid","age","Sex", "ATAC_Lane", "Fasting", "AJ_ancestry_.", "FRiP_Score", "Non_CD4_T_Proportion", "ATAC_Batch",colnames(PCS_1)[grep('PC', colnames(PCS_1))])]
covariates_PCs <- PCS_1[, !colnames(PCS_1) %in% c("WGS_sampleID", "Picard_PCT_PF_READS", "Picard_PCT_PF_READS_ALIGNED")]
#make sure there is no white space that might be treating batches as different entities 
covariates_PCs$ATAC_Batch<-gsub(" ", "", covariates_PCs$ATAC_Batch)
str(unique(covariates_PCs$ATAC_Batch)) #34 batches 

covariates_PCs <- covariates_PCs %>%
  mutate(ATAC_batch_lane = paste(ATAC_Batch, ATAC_Lane, sep="_")) %>%
  mutate(ATAC_batch_lane = as.numeric(as.factor(ATAC_batch_lane)))


# Replace 'PC' with 'ATAC_PC' in each column name that starts with 'PC'
colnames(covariates_PCs) <- gsub("^PC", "ATAC_PC", colnames(covariates_PCs))
# List of GEX_PC columns that are currently character
gex_columns <- paste("ATAC_PC", 1:13, sep="")

# Convert these character columns to numeric
covariates_PCs[gex_columns] <- lapply(covariates_PCs[gex_columns], function(x) as.numeric(as.character(x)))

# Check for any conversion problems, such as NAs introduced by coercion
sum(is.na(covariates_PCs[gex_columns]))


#----------------Merge Covariate files-----------------------

#covariates_PCs$sampleid<- gsub("g", "", covariates_PCs$sampleid)
str(unique(covariates_PCs$sampleid)) #362 samples 
covariates_PCs <- covariates_PCs[match(genotype$sampleid, covariates_PCs$sampleid), ]

covariates_PCs <- merge(covariates_PCs, genotype, by="sampleid")
dim(covariates_PCs)
#--------------Change column names and order --------------------
colnames(covariates_PCs)[colnames(covariates_PCs) == "age"] <- "Age"

#convert ATAC_lane to categorical variable
covariates_PCs$ATAC_Lane<-as.character(covariates_PCs$ATAC_Lane)
covariates_PCs$sampleid<-as.character(covariates_PCs$sampleid)
class(covariates_PCs$sampleid)
rownames(covariates_PCs) <- covariates_PCs$sampleid
covariates_PCs$sampleid <- NULL

covariates_PCs <- covariates_PCs %>% select('Age','Sex','ATAC_Lane','Fasting','AJ_ancestry_.','FRiP_Score','ATAC_batch_lane','Non_CD4_T_Proportion','ATAC_Batch','p_CD4_Naive','p_CD4_Central_Memory','p_CD4_Effector_Memory',
                                            'p_CD4_Cytotoxic','p_CD4_Treg','Weighted_Mean_GC_Content','Mean_Base_Quality','Mean_N_Content','Total_Aligned_Filtered_Reads','Picard_PF_ALIGNED_BASES','Picard_PF_HQ_ALIGNED_Q20_BASES','Picard_PF_MISMATCH_RATE',
                                            'Picard_PF_HQ_ERROR_RATE','Picard_MEAN_READ_LENGTH','Picard_MEDIAN_READ_LENGTH','Picard_STRAND_BALANCE','Picard_PCT_READS_ALIGNED_IN_PAIRS','ATAC_PC1','ATAC_PC2','ATAC_PC3','ATAC_PC4',
                                            'ATAC_PC5','ATAC_PC6','ATAC_PC7','ATAC_PC8','ATAC_PC9','ATAC_PC10','ATAC_PC11','ATAC_PC12','ATAC_PC13','ATAC_PC14','ATAC_PC15','ATAC_PC16','ATAC_PC17','ATAC_PC18','ATAC_PC19','ATAC_PC20','ATAC_PC21',
                                            'ATAC_PC22','ATAC_PC23','Genotype_PC1','Genotype_PC2','Genotype_PC3','Genotype_PC4','Genotype_PC5','Genotype_PC6','Genotype_PC7','Genotype_PC8','Genotype_PC9','Genotype_PC10',
                                            'Genotype_PC11','Genotype_PC12','Genotype_PC13','Genotype_PC14','Genotype_PC15')  # replace with your actual column names                                      
                                      
                                      
#######################################                                      
##--------------OPTIONAL --------------
##-----------Plot visualization
########################################
                                      
covariates_PCs$Sex <- ifelse(covariates_PCs$Sex == "Female", 0, 1)
covariates_PCs$ATAC_Batch <- factor(covariates_PCs$ATAC_Batch)
covariates_PCs$ATAC_Batch <- as.numeric(covariates_PCs$ATAC_Batch)

     
covariates_PCs$ATAC_Lane <- factor(covariates_PCs$ATAC_Lane)
covariates_PCs$ATAC_Lane <- as.numeric(covariates_PCs$ATAC_Lane)
                                      
covariates_PCs$ATAC_batch_lane <- factor(covariates_PCs$ATAC_batch_lane)    
covariates_PCs$ATAC_batch_lane <- as.numeric(as.factor(covariates_PCs$ATAC_batch_lane))   
                                      
covariates_PCs$Fasting <- factor(covariates_PCs$Fasting)
covariates_PCs$AJ_ancestry_. <- as.numeric(covariates_PCs$AJ_ancestry_.)
covariates_PCs$FRiP_Score  <- as.numeric(covariates_PCs$FRiP_Score )
covariates_PCs$Non_CD4_T_Proportion <- as.numeric(covariates_PCs$Non_CD4_T_Proportion)
covariates_PCs$Fasting <- as.numeric(covariates_PCs$Fasting)
        
M<-cor(covariates_PCs)
head(round(M,2))

# Save the plot to a PNG file
png( paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_cov_PC_corrplot.png') ,width = 30, height = 30, units ='cm', res=300, pointsize = 8)
corrplot(M, type="upper", tl.col="black")
dev.off()


#######################
##----SAVING TABLES----
##
## Please note: Based on visual inspection, 
## I decided that I will include Genotype and phenotype PCs
## and Sex and Age (because it looks like PCs capture most variation)
########################
                                                                            
write.table(covariates_PCs, paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_covs_PCs_all_covariates.txt'), sep = '\t', quote = F, row.names = T,
           col.names = T)

# Only PCs
write.table(covariates_PCs[,27:64], paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_covs_PCs_only.txt'), sep = '\t', quote = F, row.names = T, col.names = T) #pcs only
write.table(covariates_PCs[,c(1:2,27:64)], paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_covs_PCs_age_sex.txt'), sep = '\t', quote = F, row.names = T, col.names = T) #pcs + age and sex

#---Note: due to the fact that we have 34 batches of RNAseq, converting the "categorical" variable to numerical might be confused by the tensorQTL as to be a continous variable 1-34
#---- therefore I will use One-hot encoding: For each level of a categorical variable, one-hot encoding creates a new binary column
# Example data
batch_lane <- covariates_PCs["ATAC_batch_lane"]



# Create one-hot encoded columns
# The '-1' in the formula removes the intercept column, which is typically the first level
encoded_batches <- model.matrix(~ ATAC_batch_lane - 1, data=batch_lane)
encoded_batches <- data.frame(encoded_batches)

# View the encoded data
head(encoded_batches)

covariates_PCs_encoded <- cbind(covariates_PCs, encoded_batches)
covariates_PCs_encoded <- covariates_PCs_encoded[, !colnames(covariates_PCs_encoded) %in% "ATAC_batch_lane"]



write.table(covariates_PCs_encoded, paste0(out_dir, gsub(".txt", "", filtered_table), '_atac_covs_PCs_encoded_batch.txt'), sep = '\t', quote = F, row.names = T,
            col.names = T)



print('Finished!')


sessionInfo()