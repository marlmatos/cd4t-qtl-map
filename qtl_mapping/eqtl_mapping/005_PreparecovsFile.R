#-------------------------------------------------------
# Create covariates file for tensorQTL
# Author: Winona Oliveros
# Adpated by: Marliette Matos
#------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)


PCS_1 <- read.csv2("~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_elbow_pca_12.SCT.csv", sep = ",")
genotype <- read.table('/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v3_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.367.AF1.QC.BA.king2.hwe.ld.eigenvec', header = F) ### genotype PC are the same for all celltypes

## process and merge PC cov file
genotype <- genotype %>% 
  mutate(WGS_sampleID = paste0("g", V2 ))

genotype <- genotype[,c("WGS_sampleID","V3","V4","V5","V6","V7",  "V8",  "V9",  "V10", "V11", "V12", "V13")]


covariates_PCs <- PCS_1[,c("sampleid","age","Sex", "scRNA_batch_Dave",colnames(PCS_1)[grep('PC', colnames(PCS_1))])]

covariates_PCs <- merge(covariates_PCs, genotype, by.y="WGS_sampleID", by.x="sampleid")

#### fixing metadata for 2 pairs of samples that were swapped

covariates_PCs$Sex <- ifelse(covariates_PCs$sampleid=="g50463", "Male", covariates_PCs$Sex)
covariates_PCs$Sex <- ifelse(covariates_PCs$sampleid=="g20163", "Female", covariates_PCs$Sex)
covariates_PCs$Sex <- ifelse(covariates_PCs$sampleid=="g10045", "Female", covariates_PCs$Sex)


### save files 
# Assuming your data frame is named df and the column you want to modify is named column_name
covariates_PCs$sampleid<- gsub("g", "", covariates_PCs$sampleid)

covariates_PCs$sampleid<-as.character(covariates_PCs$sampleid)
class(covariates_PCs$sampleid)
rownames(covariates_PCs) <- covariates_PCs$sampleid
covariates_PCs$sampleid <- NULL

write.table(covariates_PCs, paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/005_PreparecovsFile/allcells_covs_PC.txt'), sep = '\t', quote = F, row.names = T,
            col.names = T)



