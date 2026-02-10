.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
options(bitmapType = "cairo")
library(ggplot2)
library(data.table)
library(tidyverse)
library(forcats)
library(dplyr)
library(cowplot)
setwd("~/cd4_qtl_paper_figures/figure_6")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

data.dir="/gcgl/sghatan/marlis_pj/coloc/"
# folder names
file_names = unlist(fread(paste0(data.dir, "preprocessed_folder_names.txt"), header = F))

count_cs <- function(file_name){
  region <- readRDS(file_name)
  if (isFALSE(region$converged)) 1L else length(region$sets$cs_index)
}

make_table_union <- function(gwas_id) {
  message("Processing: ", gwas_id)
  gwas_dir <- "/gcgnl/finemapping_autoimmune/"
  
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  total_cs <- if (length(regions) == 0) NA_integer_ else sum(sapply(regions, count_cs))
  
  eqtl_path  <- file.path(data.dir, "coloc_results/eqtl_gwas_coloc", paste0(gwas_id, "_coloc_results.csv"))
  caqtl_path <- file.path(data.dir, "coloc_results/ca_gwas_coloc", paste0("CD4T_chromatin_", gwas_id, "_coloc_results.csv"))
  
  eqtl_exists  <- file.exists(eqtl_path)
  caqtl_exists <- file.exists(caqtl_path)
  
  read_eqtl <- function(p) {
    if (!file.exists(p)) return(data.table())
    dt <- fread(p, showProgress = FALSE)
    if (nrow(dt) == 0) return(data.table())
    dt <- as.data.table(dt)
    dt[, region_cs := paste(region, idx1, sep = "_")]
    dt <- dt[pval < 1e-5 & pval_nominal < 1e-3]
    if (nrow(dt) == 0) return(data.table())
    dt <- dt[order(-PP.H4.abf)]
    dt <- dt[, .SD[1], by = .(eQTL_variant_GRC37, gene)]
    dt
  }
  
  read_caqtl <- function(p) {
    if (!file.exists(p)) return(data.table())
    dt <- fread(p, showProgress = FALSE)
    if (nrow(dt) == 0) return(data.table())
    dt <- as.data.table(dt)
    dt[, region_cs := paste(region, idx1, sep = "_")]
    dt <- dt[pval < 1e-5 & pval_nominal < 1e-3]
    if (nrow(dt) == 0) return(data.table())
    dt <- dt[order(-PP.H4.abf)]
    dt <- dt[, .SD[1], by = .(eQTL_variant_GRC37, peak)]
    dt
  }
  
  eqtl_coloc_result  <- read_eqtl(eqtl_path)
  caqtl_coloc_result <- read_caqtl(caqtl_path)
  
  # ---- define region_cs classes (set logic) ----
  eq_regions <- unique(eqtl_coloc_result$region_cs)
  ca_regions <- unique(caqtl_coloc_result$region_cs)
  
  both_regions      <- intersect(eq_regions, ca_regions)
  only_eq_regions   <- setdiff(eq_regions, ca_regions)
  only_ca_regions   <- setdiff(ca_regions, eq_regions)
  
  # ---- summary counts (region_cs counts) ----
  summary_table <- data.frame(
    trait = gwas_id,
    total_cond_indep = total_cs,
    just_eqtl  = length(only_eq_regions),
    just_caqtl = length(only_ca_regions),
    both       = length(both_regions),
    has_eqtl_file  = eqtl_exists,
    has_caqtl_file = caqtl_exists
  )
  
  # ---- build one integrated "intersections" table ----
  # annotate each ROW by the region_cs-level class, keep row detail
  eqtl_rows <- copy(eqtl_coloc_result)
  if (nrow(eqtl_rows) > 0) {
    eqtl_rows[, `:=`(
      trait = gwas_id,
      qtl_type = "eQTL",
      coloc_class = fifelse(region_cs %in% both_regions, "both", "eQTL_only")
    )]
  }
  
  caqtl_rows <- copy(caqtl_coloc_result)
  if (nrow(caqtl_rows) > 0) {
    caqtl_rows[, `:=`(
      trait = gwas_id,
      qtl_type = "caQTL",
      coloc_class = fifelse(region_cs %in% both_regions, "both", "caQTL_only")
    )]
  }
  
  intersections <- rbindlist(list(eqtl_rows, caqtl_rows), fill = TRUE)
  
  # Optional: keep only the columns you want in the supplement
  # (edit this list to taste)
  keep_cols <- c("trait", "qtl_type", "coloc_class", "region_cs",
                 "region", "idx1", "idx2", "region_GRC38", "variant_id_GRC38",
                 "GWAS_variant_GRC37", "eQTL_variant_GRC37",
                 "gene", "peak", "PP.H4.abf", "pval", "pval_nominal")
  keep_cols <- intersect(keep_cols, names(intersections))
  intersections <- intersections[, ..keep_cols]
  
  list(
    summary = summary_table,
    intersections = intersections
  )
}

combined_table <- lapply(file_names, make_table_union)

results2 <- bind_rows(lapply(combined_table, `[[`, "summary"))
fwrite(results2, "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_v2.csv")


all_intersections <- bind_rows(lapply(combined_table, `[[`, "intersections"))

all_intersections<-fread("~/cd4_qtl_paper_figures/figure_6/supplements//CD4T_all_coloc_intersections_eqtl_caqtl_both.csv")

#read peakset coordinates
peak_count <- fread("/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv") %>%
  select(c("V1", "Chr", "Start" , "End")) %>%
  mutate(peak_coordinate=paste0(Chr, ":" ,Start, "-" ,End)) %>%
  select(c("V1", "peak_coordinate"))

all_intersections <- all_intersections %>%
  left_join(
    peak_count %>% dplyr::rename(peak = V1),  # rename once, simplest
    by = "peak"
  )


fwrite(all_intersections, "~/cd4_qtl_paper_figures/figure_6/supplements//CD4T_all_coloc_intersections_eqtl_caqtl_both.csv")
