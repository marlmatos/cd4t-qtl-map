rm(list = ls())
library(data.table)
library(tidyverse)

setwd("/gchm_lab/collab/marlis_pj/coloc/")

setDTthreads(8)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("chr number")
} else {
  chr = args[1]
  name = args[2]
}

# name = "All_CD4T_cells"
# chr = 3

# Read lbf file
grc38_lbf = fread(paste0("SuSiE_finemap_lbf/",name,"/",name,"_chr",chr,"_lbfs.txt"))

# Ensure chr has "chr" prefix
grc38_lbf[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]

# Add BED-required start and end columns
grc38_lbf[, `:=`(
  start = pos - 1,
  end = pos
)]

# Reorder columns: BED3 first, then everything else
bed_formatted <- grc38_lbf[, c("chr", "start", "end", setdiff(names(grc38_lbf), c("chr", "start", "end"))), with = FALSE]

bed_file = paste0("SuSiE_finemap_lbf/",name,"/", name, "_chr", chr, ".bed")

bed_formatted[, start := format(start, scientific = FALSE)]
bed_formatted[, start := as.integer(start)]
bed_formatted[, end := format(end, scientific = FALSE)]
bed_formatted[, end := as.integer(end)]

# Write to file (no quotes, tab-delimited)
fwrite(bed_formatted, file = bed_file, sep = "\t", quote = FALSE)

print("Running bedtools intersect")

# Perform bedtools intersect
system(paste0("tail -n +2 ", bed_file, " | /nfs/sw/bedtools/bedtools-2.31.1/bin/intersectBed -wa -wb -a stdin -b /gchm_lab/collab/genome/dbsnp151_GRCh38p7/bed_chr_",
              chr,".bed.gz > bedfiles/grc38_rsids_chr",chr,".txt"))

print("Intersect done")

# read in intersection output
# Load grc38 susie lbfs with rsids
grc_38 = fread(paste0("bedfiles/grc38_rsids_chr",chr,".txt"))
colnames(grc_38) = c(colnames(bed_formatted), "chr.x", "start.grc38", "end.grc38", "rsid", "score", "strand")

# Load grc37 positions with rsids
grc_37 = fread(paste0("/gchm_lab/collab/genome/dbsnp151_GRCh37p13/bed_chr_",chr,".bed.gz"), skip = 1)
colnames(grc_37) = c("chr", "start.grc37", "end.grc37", "rsid", "score", "strand")

# Join files based on rsids
result <- grc_37[grc_38, on = "rsid"]

if("gene" %in% colnames(result)){
  result2 = result %>% 
    select("chr", "pos_GRC37" = end.grc37, "rsid", "posGRC38" = end.grc38, "variant_id_GRC38" = variant_id.x,
           "lbf_1", "lbf_2", "lbf_3", "lbf_4", "lbf_5",
           "lbf_6", "lbf_7", "lbf_8", "lbf_9", "lbf_10",
           "PIP", region_GRC38 =  region, "gene",
           "start_distance", "af", "ma_samples", "ma_count",
           "pval_nominal", "beta", "slope_se", "a0", "a1")
} else {
  result2 = result %>% 
    select("chr", "pos_GRC37" = end.grc37, "rsid", "posGRC38" = end.grc37, 
           "variant_id_GRC38" = variant_id.x,
           "lbf_1", "lbf_2", "lbf_3", "lbf_4", "lbf_5",
           "lbf_6", "lbf_7", "lbf_8", "lbf_9", "lbf_10",
           "PIP", "region", "peak",
           "start_distance", "af", "ma_samples", "ma_count",
           "pval_nominal", "beta", "slope_se", "a0", "a1")
}

fwrite(result2, paste0("SuSiE_finemap_lbf_GRC37/",name,"/",name,"_chr",chr,"_lbfs_GRC37.txt"), quote = F, row.names = F, sep = "\t")
