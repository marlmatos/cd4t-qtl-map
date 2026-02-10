##############################################################################
## Identify PCs in all cells
#
# Description: This R script was written to run a job for PCA anlaysis
##############################################################################

print(Sys.time())
start=Sys.time()

args = commandArgs(trailingOnly=TRUE)

#celltype <- args[1]
#celltype=as.numeric(ct_name)

# opt <- args[2]
# opt = as.numeric(opt)

# Import libraries
library(PCAForQTL)
library(data.table)
library(RNOmni)
library(dplyr)
library(lubridate)
library(chron)
# library(moments)
source("~/cd4_QTL_analysis/02_Gene_expression//NormalizePseudobulk.R")

set.seed(2022)

# Get expression file --------------------------------
expr=fread(paste0("~/cd4_QTL_analysis/02_Gene_expression/analysis/002_GetPseudobulk_whole/all_cells_mean_mx.txt"),header=T)
expr=as.data.frame(expr)

dim(expr)
gene=read.table(paste0("~/cd4_QTL_analysis/02_Gene_expression/analysis/002_GetPseudobulk_whole/all_cells_gene_list.txt"),header=T)
expr=as.data.frame(t(expr))
colnames(expr)=gene$gene
dim(expr)

# Remove individuals with NA
expr <- mutate_all(expr, function(x) as.numeric(as.character(x)))
na_idx=which(is.na(rowMeans((expr))))
if(length(na_idx) > 0){
  expr = expr[-na_idx,]
}

# Remove genes with NA
na_idx=which((colSums((expr))==0))
if(length(na_idx) > 0){
  expr = expr[,-na_idx]
}

dim(expr)

expr[1:5,1:5]
rnames=rownames(expr)
cnames=colnames(expr)

# Get covariate files ---------------------------------
#note: it's better to use the metadata from the single cell, here I manually select out samples that are duplicated 
#from metadata nad that have been previously removed from scRNAseq
metadata <- read.csv('~/metadata_harmonizing/results/post_WGS_QC/version_080923/revision_103023/cd4_scRNAseq_meta.v3.csv')

metadata <- filter(metadata, !scRNA_Sample_Name %in% c("A3402", "A0103", "A2004"))
# Check all the individuals exist in both files
metadata <- metadata %>% 
  mutate(WGS_sampleID.qtl = paste0("g", WGS_sampleID))

table(rnames %in% metadata$WGS_sampleID.qtl)

#Calculate Age ----------------------------------

##yeas before 1950 are in the wrong centory
metadata$DOB_formatted<-mdy(metadata$DOB)
#some of the DOBs when formatted  to YYYY from YY had assigned the wrong century
as.Date(chron(format(as.Date(metadata$DOB_formatted, "%m/%d/%y"), "%m/%d/%y"))) 
options(chron.year.expand =   #expanding the cutoff year
          function (y, cut.off = 12, century = c(1900, 2000), ...) {
            chron:::year.expand(y, cut.off = cut.off, century = century, ...)
          })
metadata$DOB_formatted <- as.Date(chron(format(as.Date(metadata$DOB_formatted, "%m/%d/%y"), "%m/%d/%y"))) #re-formatting the years with the right century
   
## Format the date collected to date format
metadata$Date_formatted <- mdy(metadata$Date) 


## calculate the age of the individual at the time of collection
metadata <-metadata %>%
  mutate(
    age = year(Date_formatted) - year(DOB_formatted))


### Read Ethnicity ------------------------
ashkenazi_ancestry <- read.csv("~/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/scripts/AJ_ancestry.csv")
dim(ashkenazi_ancestry)

colnames(ashkenazi_ancestry) <- c("WGS_sampleID", "Parent_Directory", "AJ_ancestry_%" )

metadata$WGS_sampleID <- as.character(metadata$WGS_sampleID)
ashkenazi_ancestry$WGS_sampleID <- as.character(ashkenazi_ancestry$WGS_sampleID)
metadata<- left_join(metadata, ashkenazi_ancestry, by="WGS_sampleID")
dim(metadata)

# summarize covs
metadata <- metadata[,c("WGS_sampleID.qtl", "age", "Sex", "scRNA_batch_Dave", "AJ_ancestry_%")] %>% distinct()
dim(metadata)

metadata <- metadata %>% filter(metadata$WGS_sampleID.qtl %in% rnames)

duplicates <- metadata$WGS_sampleID.qtl[duplicated(metadata$WGS_sampleID.qtl)]
print('The duplicated samples are:')                   
duplicates


# ***IMPORTANT*** Match the order of individuals
metadata = metadata[match(rnames, metadata$WGS_sampleID.qtl),]
                   
# Keep only
metadata<- metadata[!duplicated(metadata$WGS_sampleID.qtl), ]

# Set row names to WGS_sampleID                                      
rownames(metadata) <- metadata$WGS_sampleID.qtl

# Transform and scale the matrix
expr2 <- pseudobulk_scaling(expr=expr,pi0=0.9,log1p=T,scale=T,RINT=F,HVG=F)

print("Start the PCA estimation...")
print(Sys.time())
tick1 = Sys.time()

### check finite values ####
expr2[1:5,1:5]
dim(expr2)
table(is.na(expr2))
na_idx=which(is.na(colMeans((expr2))))
if(length(na_idx) > 0){
  expr2 = expr2[,-na_idx]
}
dim(expr2)

prcompResult <- prcomp(expr2, center = F, scale = F) #This should take less than a minute.
PCs <- prcompResult$x

                   
                   
                   
print(Sys.time())
tick2 = Sys.time()
z = as.numeric(tick2-tick1,units="mins")
print(paste0("Running time for PCA estimation is ",z," mins"))

### select PC based on different methods ###
## elbow ##
print(paste0("Selecting PCs with elbow method"))
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(paste0("PCs selected: ", resultRunElbow))

## BE ##
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
print(paste0("Selecting PCs with BE method"))
resultRunBE<-PCAForQTL::runBE(expr2,B=20,alpha=0.05)
print(paste0("PCs selected: ", resultRunBE$numOfPCsChosen))

## plot and compare ##
K_elbow<-resultRunElbow #12.
K_BE<-resultRunBE$numOfPCsChosen #29.
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                         titleText="all_cells")
ggplot2::ggsave(paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_PCs.SCT.pdf'),width=16,height=11,unit="cm")

# Add known covariates back
# elbow
PC_num = K_elbow
factors_df = cbind(sampleid = rownames(PCs), metadata, PCs[,1:PC_num])

#colnames(factors_df) = c(c("sampleid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age"), paste0("PCA", 1:PC_num))

write.table(factors_df, paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_elbow_pca_', PC_num, ".SCT.txt"), col.names = T, row.names = F, quote = F)
write.csv(factors_df, file = paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_elbow_pca_', PC_num, '.SCT.csv'), row.names = FALSE)

# BE
PC_num = resultRunBE$numOfPCsChosen
factors_df = cbind(sampleid = rownames(PCs), metadata, PCs[,1:PC_num])

#colnames(factors_df) = c(c("sampleid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age"), paste0("PCA", 1:PC_num))

write.table(factors_df, paste0('~/cd4_QTL_analysis/02_Gene_expression/analysis/003_Identify_PCs/allcells_pseudo_BE_pca', PC_num, ".SCT.txt"), col.names = T, row.names = F, quote = F)


print(Sys.time())
end = Sys.time()
z = as.numeric(end - start, units = "mins")
print(paste0("Running time for whole script is ", z, " mins"))
