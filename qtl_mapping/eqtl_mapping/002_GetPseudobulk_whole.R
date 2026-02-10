#################################################################
#### Script to get pseudobulk #####
### Step1: Read object
### Step2: Make the identity of the object the sample IDs from WGS
### Step3: Calculate pseudobulk averages per gene per sample using the SCTransformed counts
### Step4: Write gene list in a table
### Step5: Write sample averages in a table
##############################################
# Import packages
library(Seurat)
library(dplyr)
library(data.table)
options(future.globals.maxSize = 1e16)

### Step1: Read object
# gene expression level per scRNA_Sample_Name was calculated as the intra-scRNA_Sample_Name mean counts across cells
data.dir <- '/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/'

expression_filename <- paste0(data.dir, "001_preparing_seurat_obj_QTL/2024-08-12_all_cd4T_filt_scaled_seurat_qtl.rds")

expr <- readRDS(expression_filename)

## Normalize cells 
expr <- SCTransform(expr, vars.to.regress = "percent.mt", verbose = FALSE)

md <- expr@meta.data

## Number of cells and WGS_sampleID
print(paste("This data set includes",dim(md)[1],"all cells from ",sum(table(md$WGS_sampleID)>0)," WGS_sampleID"))

print("Number of cells per WGS_sampleID: ")
print(summary(as.numeric(table(as.character(md$WGS_sampleID)))))

# Extract a list of WGS_sampleID and their number of cells
a=unlist(table(md$WGS_sampleID))


a=as.data.frame(a)
colnames(a)<-c('WGS_sampleID', 'nCells')

#get the average number of reads per sample
average_nCounts_per_sample <- round(aggregate(nCount_RNA ~ WGS_sampleID, data = md, FUN = mean))

stats<-merge(a, average_nCounts_per_sample, by='WGS_sampleID')

write.table(stats ,paste0(data.dir, '002_GetPseudobulk_whole/all_cells_nCells_nReads.txt'),row.names=F,col.names=T,quote=F)



### Step2: Make the identity of the object the sample IDs from WGS
Idents(object=expr ) <- "WGS_sampleID"
length(levels(x=expr))


 ### Step3: Calculate pseudobulk averages per gene per sample using the SCTransformed counts
 ## Main function to calculate the pseudobulk mean expression

person.averages <- AverageExpression(object= expr, assays = "SCT", features = NULL,
                                   return.seurat = FALSE, group.by = "WGS_sampleID", layer = "counts",
                                    use.scale = FALSE, use.counts = FALSE, verbose = TRUE)

person.averages=as.data.frame(person.averages[[1]])
 
gene=rownames(person.averages)

# ## Save results
# ### Step4: Write gene list in a table
write.table(as.data.frame(gene),paste0(data.dir, '002_GetPseudobulk_whole/all_cells_gene_list.txt'),row.names=F,col.names=T,quote=F)
# 
# 
# ### Step5: Write sample averages in a table
write.table(person.averages,paste0(data.dir, '002_GetPseudobulk_whole/all_cells_mean_mx.txt'),row.names=F,col.names=T,quote=F)

sessionInfo()
