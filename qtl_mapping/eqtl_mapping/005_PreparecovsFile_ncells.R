#-------------------------------------------------------
# Create covariates file for tensorQTL
# Author: Winona Oliveros
# Adpated by: Marliette Matos
#------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)
library(Hmisc)
library(corrplot)
library(ggplot2)
library(RColorBrewer)


PCS_1 <- read.csv2("~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_elbow_pca_12.SCT.csv", sep = ",")
genotype <- read.table('/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.ld.eigenvec', header = F) ### genotype PC are the same for all celltypes

## process and merge PC cov file
genotype <- genotype %>% mutate(WGS_sampleID = paste0("g", V2 ))
colnames(genotype)[colnames(genotype) == "V2"] <- "WGS_sampleID"
genotype <- genotype[, !colnames(genotype) %in% "V1"]
#-----------subset to adequate number of genotype PCS
genotype <- genotype[,c("WGS_sampleID","V3","V4","V5","V6","V7","V8", "V9", "V10","V11","V12","V13")]

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


#-----------read other source of variation such as ncells and nRNA_count
ncells<-read.table("~/cd4_QTL_analysis/02_Gene_expression/analysis/002_GetPseudobulk_whole/all_cells_nCells_nReads.txt", header = T)


covariates_PCs <- PCS_1[,c("sampleid","age","Sex", "scRNA_batch_Dave",colnames(PCS_1)[grep('PC', colnames(PCS_1))])]
# Replace 'PC' with 'GEX_PC' in each column name that starts with 'PC'
colnames(covariates_PCs) <- gsub("^PC", "GEX_PC", colnames(covariates_PCs))
# List of GEX_PC columns that are currently character
gex_columns <- paste("GEX_PC", 1:12, sep="")

# Convert these character columns to numeric
covariates_PCs[gex_columns] <- lapply(covariates_PCs[gex_columns], function(x) as.numeric(as.character(x)))

# Check for any conversion problems, such as NAs introduced by coercion
sum(is.na(covariates_PCs[gex_columns]))

#make sure there is no white space that might be treating batches as different entities 
covariates_PCs$scRNA_batch_Dave<-gsub(" ", "", covariates_PCs$scRNA_batch_Dave)
str(unique(covariates_PCs$scRNA_batch_Dave)) #34 batches 


#----------------Merge Covariate files-----------------------
covariates_PCs$sampleid<- gsub("g", "", covariates_PCs$sampleid)
str(unique(covariates_PCs$sampleid)) #364 batches 

covariates_PCs <- merge(covariates_PCs, genotype, by.y="WGS_sampleID", by.x="sampleid")

covariates_PCs<-merge(covariates_PCs, ncells, by.x="sampleid", by.y="WGS_sampleID" )

#--------------Change column names and order --------------------
colnames(covariates_PCs)[colnames(covariates_PCs) == "scRNA_batch_Dave"] <- "scRNAseq_batch"
colnames(covariates_PCs)[colnames(covariates_PCs) == "age"] <- "Age"

#-------------reorganize columns -------------------------------
total_cols <- ncol(covariates_PCs)
new_order <- c(1:4, (total_cols-1):total_cols, 5:(total_cols-2))
covariates_PCs <- covariates_PCs[, new_order]


covariates_PCs$sampleid<-as.character(covariates_PCs$sampleid)
class(covariates_PCs$sampleid)
rownames(covariates_PCs) <- covariates_PCs$sampleid
covariates_PCs$sampleid <- NULL

write.table(covariates_PCs, paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/005_PreparecovsFile_ncells/allcells_covs_PC_nFeatures_ncells.txt'), sep = '\t', quote = F, row.names = T,
           col.names = T)


#---Note: due to the fact that we have 34 batches of RNAseq, converting the "categorical" variable to numerical might be confused by the tensorQTL as to be a continous variable 1-34
#---- therefore I will use One-hot encoding: For each level of a categorical variable, one-hot encoding creates a new binary column
# Example data
batch <- covariates_PCs["scRNAseq_batch"]


# Create one-hot encoded columns
# The '-1' in the formula removes the intercept column, which is typically the first level
encoded_batches <- model.matrix(~ scRNAseq_batch - 1, data=batch)
encoded_batches <- data.frame(encoded_batches)

# View the encoded data
head(encoded_batches)

covariates_PCs_encoded <- cbind(covariates_PCs, encoded_batches)
covariates_PCs_encoded <- covariates_PCs_encoded[, !colnames(covariates_PCs_encoded) %in% "scRNAseq_batch"]

write.table(covariates_PCs_encoded, paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/005_PreparecovsFile_ncells/allcells_covs_PC_nFeatures_ncells_encoded_batch.txt'), sep = '\t', quote = F, row.names = T,
            col.names = T)


##--------------OPTIONAL --------------
##-----------Plot visualization
covariates_PCs$Sex <- ifelse(covariates_PCs$Sex == "Female", 0, 1)
covariates_PCs$scRNAseq_batch <- factor(covariates_PCs$scRNAseq_batch)
covariates_PCs$scRNAseq_batch <- as.numeric(covariates_PCs$scRNAseq_batch)

M<-cor(covariates_PCs)
head(round(M,2))

# Save the plot to a PNG file
png("~/cd4_QTL_analysis/02_Gene_expression/analysis/005_PreparecovsFile_ncells/allcells_cov_PCs.png",width = 15, height = 15, units ='cm', res=300, pointsize = 8)
corrplot(M, type="upper", tl.col="black")
dev.off()





print('Finished!')


