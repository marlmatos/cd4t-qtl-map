########################################################
# Process scRNA-seq data for tensorQTL
# Author: Winona Oliveros
# Adapted by: Marliette Matos
################################################

#library(here)
library(data.table)
library(dplyr)
library(tidyr)

# Arguments ------------
annot_file <- '/gchm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf'
out_file <- "/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/004_ProcessExpression/allcells_pseudo_cells_mean_mx.bed"

# tpm_threshold <- 0.1
# count_threshold <- 6
# sample_frac_threshold <- 0.2

## format
# chr, start, end, phenotype_id, with the remaining columns corresponding to samples

# Read in data -----------

# Get expression file --------------------------------
expr=fread(paste0("/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/002_GetPseudobulk_whole/all_cells_mean_mx.txt"),header=T)
expr=as.data.frame(expr)

dim(expr)
gene=read.table(paste0("/gchm/cd4_QTL_analysis/02_Gene_expression/analysis/002_GetPseudobulk_whole/all_cells_gene_list.txt"),header=T)
expr=as.data.frame(t(expr))
colnames(expr)=gene$gene

# Dropping the "g" in the sample names obtained during the psedobulk
rownames(expr) <- gsub("g", "", rownames(expr))

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

expr <- t(expr)
pi0 = rowSums(expr==0)/ncol(expr)

# Remove genes with pi0 > 0.9
expr = expr[which(pi0 <= 0.9), ]
dim(expr)

# log(x+1) and standardization
# try without
expr = apply(expr, 2, function(x)log(x+1))
expr = as.data.frame(t(scale(t(expr))))
dim(expr)
head(expr[,c(1:5)])

## subset genes
rownames(gene) <- gene$gene
gene = as.data.frame(gene[rownames(expr),])

colnames(gene) <- c('gene')
expr = cbind(gene, expr)

# Annotation file
annot <- read.table(annot_file, header = F, sep = "\t", stringsAsFactors = F)
## Keep only genes from chr1-22
annot <- annot[annot$V1 %in% c(paste0("chr", 1:22)), ]
annot <- annot[annot$V3 %in% "gene", ]

#splitting column V9 to get gene IDs
annot$gene_id <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub("gene_id ", "", unlist(strsplit(x, ";"))[1]), "[.]"))[1]
})
#splitting column V9 to get gene names
annot$gene_name <- sapply(annot$V9, function(x) {
  unlist(strsplit(sub(".*gene_name ", "", unlist(strsplit(x, ";"))[3]), "[.]"))[1]
})
## Add start (TSS -1) and end (TSS)
## Note: if strand == +, then (start - 1, start)
## Note: if strand == -, then (end -1, end)
annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
annot$end <- ifelse(annot$V7 %in% "+", annot$V4, annot$V5)

# try instead of TSS input gene coordinates
# annot$start <- ifelse(annot$V7 %in% "+", annot$V4 - 1, annot$V5 - 1)
# annot$end <- ifelse(annot$V7 %in% "+", annot$V5, annot$V4)
# sort annot file by chr and start
annot$chr_number <- as.numeric(sub("chr", "", annot$V1))
annot <- annot[order(annot$chr_number, annot$start),]


# Extract the overlapped genes
id = Reduce(intersect, list(annot$gene_name, expr$gene))
annot = annot[match(id, annot$gene_name), ]
expr = expr[match(id, expr$gene), ]

# Combine location and expression
annot$chr <- annot$V1
data = cbind(annot[,c("chr","start","end", "gene_name")], expr[,-1])
data = data[order(data$chr, data$start), ]
colnames(data)[1]="#chr"

# Remove NAs
data = na.omit(data)

write.table(data, out_file, row.names = F, col.names = T, quote = F, sep = "\t")


cat("Done!", fill = T)
