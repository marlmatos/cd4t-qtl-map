#################################################################
#### Script to convert Seurat Object to h5ad format            ##
#### Author: Marliette Matos                                   ##
#### Project: CD4 Aging: CellRegMap of CD4 T cells             ##
#################################################################

.libPaths('~/backup_R/x86_64-pc-linux-gnu-library')
#----------load libraries---------------------
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk, lib.loc="/nfs/sw/miniconda3/miniconda3-3.22.0/envs/R-4.2.3/lib/R/library"))

print('Date ran:')
Sys.time()


#----------Read Seurat Object-----------------
## Note: The following object has been previously filtered and normalized as per QC standards 
## filtered for 'singlet' status from SouporCells & RNAcount.max=10700 & 
## RNAcount.min=500, Feature_max=2600 & Feature_min=500, mitochrial transcripts < 15% 
## filtering for CD4T cells only

root_path="~/scRNAseq/002_seurat_preprocessing/analysis/003_clustering_obj/"
dir_path="~/cd4_CellRegMap/001_preprocessing/results_01312025/prep_phenotype_vector/"

#----------------Read the gene expression data----------------
obj <- readRDS(paste0(root_path, "cd4_T_only_filt_azimuth_harmony.rds"))
obj

#---------------PREP OBJ for common samples ------------------

#Read samples that will be included in WGS
wgs_samples <- read.delim("~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.annot.fam", header = F)
wgs_samples<-wgs_samples$V2

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
write.table(common_samples_tab, file = paste0(dir_path, "/common_samples_gfgex_wgs.in.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

#Subsetting obj for just the samples that are common between GEX and WGS
obj<-subset(x = obj, idents = common_samples)

print("Number of WGS after subset:")
str(as.data.frame(table(Idents(obj))))

#--------Add additional metadata regarding batch
seq_batch <- read.csv(paste0(root_path, "cell_barcodes_with_lane.csv"))
rownames(seq_batch) <- seq_batch$cell_barcodes

#get the obj metadata
metadata <- obj@meta.data
dim(metadata)

# Subset seq_batch to keep only rows that are present in the metadata row names
seq_batch <- seq_batch[rownames(metadata), , drop = FALSE]
dim(seq_batch)
# Merge metadata with the subsetted seq_batch by row names
merged_data <- cbind(metadata, scRNAseq_batch_combo = seq_batch[, 7])
dim(merged_data)
# update the Seurat object's metadata
obj@meta.data <- merged_data

######################################-----------SPLIT OBJ PER CELLLTYPE and convert to AnnData -------------------------
#remove the scale.data slot because it interacts with the conversion to anndata, and it only tranfers a counts matrix with the highly variable genes and not all genes
#note: the only way to remove scale.data slot is by making a new object for now, different versions of seurat 
#note: this method creates a warning that the data and scale.data are empty, but that is because the normalized counts in slot "data" are 
#now put in slot "counts" in the new seurat object. I manually checked that these get incorporated the right way.

print("creating new object with just the raw counts")
# Create a new Seurat object using the normalized data (not raw counts)
export <- CreateSeuratObject(counts = GetAssayData(obj, assay = "RNA", layer = "counts"), 
                             meta.data = obj@meta.data)

#Transfer UMAP and/or other dimensions, two methods here, second one with [[""]] is better
#export@reductions$pca = obj@reductions$pca
#export@reductions$harmony = obj@reductions$harmony
#export[["umap"]] = obj[["umap"]]


#Transfer variable features if needed
#VariableFeatures(export) = VariableFeatures(obj)

print("Cell labels:") 
unique(export@meta.data$cd4_subtypes_l1) #list of cells present 
export@meta.data$cd4_subtypes_l1 <- gsub(" ", "_", export@meta.data$cd4_subtypes_l1) #get rid of any gaps in the naming if present

cell_types_lv1 <- unique(export@meta.data$cd4_subtypes_l1) #list of cellstypes present 
writeLines(cell_types_lv1, paste0(dir_path, "cell_types_lv1.txt"))

#split_seurat <- SplitObject(obj, split.by = "cd4_subtypes_l1")

#split the gene expression object for each cell type and save them as h5ad

# for (cell_type in cell_types_lv1){
#    print(cell_type)
#       #seurat_obj <- split_seurat[[cell_type]]
#     Idents(object = export) <- "cd4_subtypes_l1"
#     seurat_obj <- subset(x = export, idents = cell_type)
    
#     # Export PCA embeddings
#     #write.csv(Embeddings(seurat_obj, "pca"), file = paste0(dir_path, "/", cell_type, "_pca_embeddings.csv"))

#     # Export harmony integration embeddings
#     #write.csv(Embeddings(seurat_obj, "harmony"), file = paste0(dir_path, "/", cell_type, "_harmony_embeddings.csv"))
    

#     #---------------Convert Seurat Object to AnnData format ----------------------
#     #note: this does't work if the the file to create already exists, unless overwrite =TRUE

#     # Convert and save the Seurat object to AnnData format
#     seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay") 

#     SaveH5Seurat(seurat_obj, filename = paste0(dir_path, "/", cell_type, "_cd4_aging_filt_sce.h5Seurat"), overwrite = TRUE)
#     Convert(paste0(dir_path, "/", cell_type, "_cd4_aging_filt_sce.h5Seurat"), dest = "h5ad", overwrite = TRUE)
    
#     }


# Convert and save the Seurat object to AnnData format
export[["RNA"]] <- as(object = export[["RNA"]], Class = "Assay") 
#save the gene expression object for all cells as h5ad
SaveH5Seurat(export, filename = paste0(dir_path, "allcells_cd4_aging_filt_sce.h5Seurat"), overwrite = TRUE)
Convert(paste0(dir_path, "allcells_cd4_aging_filt_sce.h5Seurat"), dest = "h5ad", overwrite = TRUE)


sessionInfo()