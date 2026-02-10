#!/usr/bin/env Rscript
.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4",
            "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

options(bitmapType = "cairo")

# ============================================================
# Inputs
# ============================================================
data.dir   <- "/gcgl/sghatan/marlis_pj/coloc/"
gwas_dir   <- "/gcgnl/finemapping_autoimmune/"
out_file   <- "~/cd4_qtl_paper_figures/figure_6/data/SuppTable_ChromBPNet_GWAS_CS_intersections_withVariantIDs.tsv.gz"

coloc_table_path <- "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_v2.csv"
folder_names_txt <- file.path(data.dir, "preprocessed_folder_names.txt")

cbpnet_status_path <- "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_QTL_status.csv"
cbpnet_hg37_path   <- "~/cd4_qtl_paper_figures/figure_6/data/liftover_cbpnet_rsID/final/cbpnet_hg38_to_hg19_rsIDjoin.tsv.gz"

if (!file.exists(coloc_table_path)) stop("Missing: ", coloc_table_path)
if (!file.exists(folder_names_txt)) stop("Missing: ", folder_names_txt)
if (!file.exists(cbpnet_status_path)) stop("Missing: ", cbpnet_status_path)
if (!file.exists(cbpnet_hg37_path)) stop("Missing: ", cbpnet_hg37_path)

# ============================================================
# GWAS filter: only those present in coloc summary
# ============================================================
coloc_table <- fread(coloc_table_path)
coloc_table[, trait := sub("_preprocessed$", "", trait)]
trait_ids <- unique(coloc_table$trait)

fn_dt <- fread(folder_names_txt, header = FALSE)
file_names <- unlist(fn_dt, use.names = FALSE)        # "..._preprocessed"
file_base  <- sub("_preprocessed$", "", file_names)   # comparable IDs

keep_idx <- file_base %in% trait_ids
file_names_filt <- file_names[keep_idx]

message("GWAS folders retained: ", length(file_names_filt))

# ============================================================
# Load ChromBPNet liftover + attach QTL status category
# ============================================================
cbpnet_status <- fread(cbpnet_status_path)  # must contain: variant_id, category
cbpnet_hg37   <- fread(cbpnet_hg37_path)

# Try to detect rsID column if present (common names)
rs_col <- intersect(names(cbpnet_hg37), c("rsid", "rsID", "RSID", "snp", "SNP"))
rs_col <- if (length(rs_col) >= 1) rs_col[1] else NA_character_

cbpnet_hg37 <- as.data.table(cbpnet_hg37) %>%
  left_join(cbpnet_status, by = "variant_id") %>%
  as.data.table()

# Keep only successfully lifted
cbpnet_hg37 <- cbpnet_hg37[!is.na(hg19_chr) & !is.na(hg19_pos)]

# Build hg19 chr:pos keys (strip "chr")
cbpnet_hg37[, chrpos := paste0(gsub("^chr", "", hg19_chr), ":", hg19_pos)]

# Keep mapping: chrpos -> cbpnet variant identity (+ category, + rsid if available)
keep_cols <- c("chrpos", "variant_id", "category")
if (!is.na(rs_col)) keep_cols <- c(keep_cols, rs_col)

cbpnet_pos_var <- unique(cbpnet_hg37[, ..keep_cols])
setDT(cbpnet_pos_var)

message("Total cbpnet (hg19 chr:pos -> variant_id) rows: ", nrow(cbpnet_pos_var))

autosomes <- as.character(1:22)

# ============================================================
# Helpers
# ============================================================
strip_to_chrpos <- function(x) sub(":[ACGT]+:[ACGT]+$", "", x)

count_cs_region <- function(file_name) {
  region <- readRDS(file_name)
  if (isFALSE(region$converged)) 1L else length(region$sets$cs_index)
}

# ============================================================
# Region → per-CS overlap with CBPNet *variant IDs*
# ============================================================
membership_region_with_ids <- function(region, gwas_id, region_id, cbpnet_pos_var) {
  
  if (isFALSE(region$converged) || is.null(region$sets$cs) || length(region$sets$cs) == 0) {
    return(data.table())
  }
  
  susie_ids <- names(region$pip)
  if (is.null(susie_ids)) stop("Region has no variant names in 'pip'.")
  
  region_chr <- unique(sub(":.*", "", susie_ids))
  region_chr <- region_chr[1]
  
  cs_list  <- region$sets$cs
  cs_names <- names(cs_list)
  if (is.null(cs_names)) cs_names <- paste0("L", seq_along(cs_list))
  
  # If not autosomal: return CS rows as ineligible
  if (!(region_chr %in% autosomes)) {
    return(data.table(
      gwas_id   = gwas_id,
      region_id = region_id,
      chr       = region_chr,
      cs_id     = cs_names,
      eligible  = FALSE,
      has_cbpnet = FALSE,
      cs_category = NA_character_,
      n_cbpnet_variants = 0L,
      cbpnet_variant_ids = NA_character_,
      cbpnet_rsids = NA_character_,
      cbpnet_categories = NA_character_,
      susie_cs_chrpos = NA_character_
    ))
  }
  
  # Build variant-level table for this region
  dt_var <- data.table(
    idx      = seq_along(susie_ids),
    susie_id = susie_ids,
    chrpos   = strip_to_chrpos(susie_ids)
  )
  
  # Join CBPNet IDs at chrpos (can be many-to-many)
  dt_join <- cbpnet_pos_var[dt_var, on = "chrpos", nomatch = 0L, allow.cartesian = TRUE]
  # dt_join has dt_var cols plus: variant_id, category, rsid?
  
  # Summarise per CS
  out <- rbindlist(lapply(seq_along(cs_list), function(k) {
    idx_vec <- cs_list[[k]]
    
    # chr:pos of all SuSiE variants in this CS (nice for auditing)
    cs_chrpos_all <- unique(dt_var[idx %in% idx_vec, chrpos])
    
    cs_rows <- dt_join[idx %in% idx_vec]
    has_cb  <- nrow(cs_rows) > 0
    
    if (!has_cb) {
      return(data.table(
        gwas_id   = gwas_id,
        region_id = region_id,
        chr       = region_chr,
        cs_id     = cs_names[k],
        eligible  = TRUE,
        has_cbpnet = FALSE,
        cs_category = NA_character_,
        n_cbpnet_variants = 0L,
        cbpnet_variant_ids = NA_character_,
        cbpnet_rsids = NA_character_,
        cbpnet_categories = NA_character_,
        susie_cs_chrpos = paste(cs_chrpos_all, collapse = ";")
      ))
    }
    
    # Unique CBPNet variants overlapping this CS
    cb_ids  <- unique(cs_rows$variant_id)
    cb_cats <- unique(cs_rows$category)
    
    # Optional rsIDs
    cb_rs <- NA_character_
    if (!is.na(rs_col) && rs_col %in% names(cs_rows)) {
      cb_rs <- paste(unique(cs_rows[[rs_col]]), collapse = ";")
    }
    
    # Precedence category for CS
    cs_cat <- if ("QTL_shared" %in% cb_cats) {
      "QTL_shared"
    } else if ("ChromBPNet_specific" %in% cb_cats) {
      "ChromBPNet_specific"
    } else {
      NA_character_
    }
    
    data.table(
      gwas_id   = gwas_id,
      region_id = region_id,
      chr       = region_chr,
      cs_id     = cs_names[k],
      eligible  = TRUE,
      has_cbpnet = TRUE,
      cs_category = cs_cat,
      n_cbpnet_variants = length(cb_ids),
      cbpnet_variant_ids = paste(cb_ids, collapse = ";"),
      cbpnet_rsids = cb_rs,
      cbpnet_categories = paste(cb_cats, collapse = ";"),
      susie_cs_chrpos = paste(cs_chrpos_all, collapse = ";")
    )
  }), use.names = TRUE, fill = TRUE)
  
  out
}

# ============================================================
# GWAS → all regions → bind per-CS table
# ============================================================
make_intersections_table <- function(gwas_folder, cbpnet_pos_var) {
  
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_folder),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  if (length(regions) == 0L) return(data.table())
  
  total_cs_all <- sum(sapply(regions, count_cs_region))
  
  cs_tbl <- rbindlist(lapply(regions, function(f) {
    region <- readRDS(f)
    region_id <- sub("\\.rds$", "", basename(f))
    membership_region_with_ids(region, gwas_folder, region_id, cbpnet_pos_var)
  }), use.names = TRUE, fill = TRUE)
  
  cs_tbl[, total_cs_all := total_cs_all]
  cs_tbl
}

# ============================================================
# Run and write
# ============================================================
all_cs_intersections <- rbindlist(
  lapply(file_names_filt, make_intersections_table, cbpnet_pos_var = cbpnet_pos_var),
  use.names = TRUE, fill = TRUE
)

message("Rows written: ", nrow(all_cs_intersections))
fwrite(all_cs_intersections, out_file, sep = "\t")
message("Wrote: ", out_file)

all_cs_intersections<-fread("cd4_qtl_paper_figures/figure_6/supplements/SuppTable_ChromBPNet_GWAS_CS_intersections_withVariantIDs.tsv")
all_cs_intersections <- all_cs_intersections %>% filter(has_cbpnet==TRUE) 

fwrite(all_cs_intersections, "cd4_qtl_paper_figures/figure_6/supplements/SuppTable_ChromBPNet_GWAS_CS_intersections_withVariantIDs.tsv", sep = "\t")
