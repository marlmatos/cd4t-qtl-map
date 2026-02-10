#################################################################
## This script is to prepare the seurat object for QTL analysis
## previous: "~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/scripts/sam_MAF1/v2_removing_outliers/002.v2_calculating_pcs_MAF1-ran.ipynb"
## step 1: read seurat object
## step 2: set the identity of the object to the Azimuth labeling of the cells
## step 3: Subset the obj to just CD4+ cells
## step 4: Read the .fam file from the pre-QC of the WGS to know get the sample IIDs of all unrelated samples
## step 5: Subset obj to contain only the common samples between both
## step 6: write a table in plink format witrh the samples iids common between the datasets
## step 7: Save Seurat obj
#################################################################

suppressPackageStartupMessages(library(Seurat))


root_dir="/gchm/scRNAseq/002_seurat_preprocessing/analysis/002_creating_seurat_objects/merged_object/"
## step 1
# Path of gene expression data out of CellRanger
obj <- readRDS(paste0(root_dir, "cd4_all_cells_filt_dedup_azimuth.rds"))

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj

## step 2:
#filter objects for just CD4 T cells with a Azimuth assignment probablity of > 0.8
Idents(object = obj) <- "predicted.celltype.l1"

## step 3:
obj <-subset(obj, idents = "CD4 T", predicted.celltype.l1.score>0.8)

## step 4:
#Read samples that will be included in WGS
wgs_samples <- read.delim("/gchm/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.annot.fam", header = F)
wgs_samples<-wgs_samples$V2

# step 5: 
#Extracting common samples between gene expression and whole genome sequencing experiments
# Set identity classes to an existing column in meta data
Idents(object = obj) <- "WGS_sampleID"

gex_samples <- as.data.frame(table(Idents(obj)))

print("The number of common samples between WGS and scGEX is:")
summary(wgs_samples %in% gex_samples$Var1)

gex_samples<-as.vector(gex_samples$Var1)

 
common_samples <-wgs_samples[wgs_samples %in% gex_samples]


# saving the the sample list as a plink compatible file
common_samples_tab<-as.data.frame(common_samples)
common_samples_tab$FIID <- "0"
colnames(common_samples_tab) <- c("IID", "FIID")
# Swap the order of columns while keeping the contents intact
common_samples_tab <- common_samples_tab[, c("FIID", "IID")]

print("Saving Plink sample IIDs")
write.table(common_samples_tab, file = "/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/001_preparing_seurat_obj_QTL/003_common_samples_gfgex_wgs.in.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

obj<-subset(x = obj, idents = common_samples)

print("Number of WGS after subset:")
str(as.data.frame(table(Idents(obj))))

# Save merged Seurat object
saveRDS(obj, paste0("/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/001_preparing_seurat_obj_QTL/",Sys.Date(), "_all_cd4T_filt_scaled_seurat_qtl.rds"))
