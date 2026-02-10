#------------------------
# Description: This script is to prepare the ATACseq peak quantification matrix for caQTL analysis
#
# Author: Marliette Matos
# Date: 01/26/2025
# 
#--------------------------

#####
##  NOTE: The genotype file used in this analysis was QCed and preprocessed in
##  "cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/scripts/sam_MAF5/v4_removing_PC_outliers" and
## it is the same used for eQTL mapping for this project
#####

# Libraries
.libPaths("/gpfs/commons/home/mmatos/R/x86_64-pc-linux-gnu-library")
library(edgeR)
library(DESeq2)
library(matrixStats)
library(qsmooth)
library(dplyr)
library(GenomicAlignments)
set.seed(2024)

# Directories
root_dir <- "/gpfs/commons/home/mmatos/scRNAseq/002_seurat_preprocessing/analysis/004_creating_seurat_objects_filtered_annotated/"
out_dir <- "/gpfs/commons/home/mmatos/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/001_peaks/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


##
##--------Read Phenotype data tables 
## Paths
peaks_counts <- "/gpfs/commons/home/mmatos/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"
peak_metadata <- "/gpfs/commons/home/mmatos/ATAC-seq_analysis/diff_accesibility_ana/results/metadata/cd4_atac_metadata_scrna_union.csv"

## Read data
peak_metadata <- read.csv(peak_metadata, row.names = 1)

## Get the gwnomic coordinates for each peak
ranges.table <- read.csv2(peaks_counts, sep=",")
ranges.table<-as.data.frame(ranges.table)
rownames(ranges.table) <- ranges.table$X
chroms_to_keep <- paste0("chr", 1:22)
ranges.table_filtered <- ranges.table %>% 
  filter(Chr %in% chroms_to_keep)

#drop the. coordinates to just get peaks
count.table<- ranges.table_filtered[,7:ncol(ranges.table_filtered),drop=FALSE]


print("Reading MetaData")
print(dim(peak_metadata))
print(peak_metadata[1:5, 1:5])

print("Reading Peak Counts")
print(dim(count.table))
print(count.table[1:5, 1:5])


###------------------------
# GC-aware smooth quantile Normalization
# get GC content
rn <- rownames(ranges.table_filtered)
sn <- ranges.table_filtered$Chr
start <- as.numeric(ranges.table_filtered$Start)
end <- as.numeric(ranges.table_filtered$End)
gr <- GRanges(seqnames=sn, ranges=IRanges(start, end), strand="*", mcols=data.frame(peakID=rn))
ff <- FaFile("/gpfs/commons/home/mmatos/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa")


peakSeqs <- getSeq(x=ff, gr)
gcContentPeaks <- letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
gcGroups <- Hmisc::cut2(gcContentPeaks, g=20)
mcols(gr)$gc <- gcContentPeaks


###----------------------------
# Filtering 
keep <- rowSums(cpm(count.table) >= 2) >= (0.1 * ncol(count.table))
count.table.keep <- count.table[keep, ]
gc <- gcContentPeaks[keep]

# GC aware starndadization
qs_norm_gc <- qsmoothGC(object = count.table.keep, gc=gc, 
                        group_factor = peak_metadata$Group)  #grouping criteria needs to specified to avoid over correcting for biological  effects, this is more specific for different celltype, doesnt make much difference by group because celltypes are the same across individuals


# Get the numeric matrix from qsmooth object
qsmooth_matrix <- qsmoothData(qs_norm_gc)


#convert normalized counts to cpm
cpmqsmooth <- edgeR::cpm(qsmooth_matrix, 
                         normalized.lib.sizes=FALSE,  #no need for normalizatrion because it has been already smooth quantile nromalized
                         log=FALSE, #log normalize 
                         prior.count=0.25)

# Read WGS samples
wgs_samples <- read.delim("/gpfs/commons/home/mmatos/cd4_QTL_analysis/01_genotype_snps_covar/02_genotype_covariates/analysis/002.v4_calculating_pcs_MAF5/CD4_all_chr_ashkenazi.364.AF1.QC.BA.king2.hwe.annot.fam", header = F)
wgs_samples <- wgs_samples$V2

# Filter samples
common_samples <- peak_metadata$ATAC_Sample_Name[peak_metadata$WGS_sampleID %in% wgs_samples]
common_samples2 <- peak_metadata$WGS_sampleID[peak_metadata$WGS_sampleID %in% wgs_samples]
atac_sample_names <- peak_metadata$ATAC_Sample_Name[peak_metadata$WGS_sampleID %in% wgs_samples]
write.table(atac_sample_names, file = paste0(out_dir, "qtl_atac_sample_names_convention.txt"), sep = "\t", quote = FALSE, col.names = NA)

# Save sample list
common_samples_tab <- data.frame(FIID = "0", IID = common_samples2)
write.table(common_samples_tab, file = paste0(out_dir, "001_common_samples_atac_wgs.in.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Subset peak counts
cpmqsmooth <- cpmqsmooth[, common_samples]
colnames(cpmqsmooth) <- peak_metadata$WGS_sampleID[match(colnames(cpmqsmooth), peak_metadata$ATAC_Sample_Name)]

print("Number of samples in peak_file after WGS subset:")
print(dim(cpmqsmooth))


#save peaks
write.table(cpmqsmooth, file = paste0(out_dir, "filtered_qsmooth_norm_cpm.txt"), sep = "\t", quote = FALSE, col.names = NA)
            
# Save peak names
write.table(rownames(cpmqsmooth), file = paste0(out_dir, "peak_names.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

print("Analysis complete. Results saved as tab-separated text files in the output directory.")

sessionInfo()