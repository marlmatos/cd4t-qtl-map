#!/usr/bin/env Rscript
.libPaths(c(
  "/gchm/R/x86_64-pc-linux-gnu-library/4.4",
  "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"
))
options(bitmapType = "cairo")

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(forcats)

setwd("~/cd4_qtl_paper_figures/figure_6")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

## ----------------------------------------------------------------------
## 0) Inputs
## ----------------------------------------------------------------------

#######read qtl-gwas colocalization summary table
coloc_table = fread("~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_v2.csv")
# Clean `trait` column (remove "_preprocessed")
coloc_table$trait <- gsub("_preprocessed$", "", coloc_table$trait)


data.dir <- "/gpfs/commons/groups/lappalainen_lab/sghatan/marlis_pj/coloc/"

## 1) Read folder names as a simple character vector
fn_dt <- fread(file.path(data.dir, "preprocessed_folder_names.txt"),
               header = FALSE)

# Flatten all columns to a single character vector
file_names <- unlist(fn_dt, use.names = FALSE) # file_names looks like: "2021_34594039_AD_EUR-EAS_preprocessed", etc.

## 2) Get the trait IDs used in coloc_table
trait_ids <- unique(coloc_table$trait)
# e.g. "2021_34594039_AD_EUR-EAS"

## 3) Strip "_preprocessed" from file_names to get comparable IDs
file_base <- sub("_preprocessed$", "", file_names)
# e.g. "2021_34594039_AD_EUR-EAS"

## 4) Keep only those GWAS that are in coloc_table$trait
keep_idx <- file_base %in% trait_ids

file_names_filt  <- file_names[keep_idx]   # folders you can pass to CBPNet
gwas_ids_filt    <- file_base[keep_idx]    # matching trait-style IDs

##sanity checks
setdiff(trait_ids, gwas_ids_filt)  # traits in coloc_table but no folder
setdiff(gwas_ids_filt, trait_ids)  # folders with no coloc summary


# ChromBPNet QTL status (has 'variant_id', 'category', etc.)
cbpnet_status <- fread(
  "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_QTL_status.csv"
)

# Example GWAS id (for quick test)
gwas_id <- file_names[101]

# ChromBPNet hg38→hg19 results with rsID join
chrompbnet_hg37 <- "~/cd4_qtl_paper_figures/figure_6/data/liftover_cbpnet_rsID/final/cbpnet_hg38_to_hg19_rsIDjoin.tsv.gz"

if (!file.exists(chrompbnet_hg37)) {
  stop("ChromBPNet hg37 result file not found: ", chrompbnet_hg37)
}

## ----------------------------------------------------------------------
## 1) Preprocess ChromBPNet variants → hg19 chr:pos + category
## ----------------------------------------------------------------------

chrompbnet_hg37_result <- fread(chrompbnet_hg37)

chrompbnet_hg37_result <- chrompbnet_hg37_result %>%
  left_join(cbpnet_status, by = "variant_id")

# Keep only successfully lifted variants
chrompbnet_hg37_result <- chrompbnet_hg37_result[
  !is.na(hg19_chr) & !is.na(hg19_pos)
]

# Build chr:pos IDs (strip "chr" prefix to match SuSiE style)
chrompbnet_hg37_result[, chrpos := paste0(
  gsub("^chr", "", hg19_chr), ":", hg19_pos
)]

# Position ↔ category mapping
# (can be multiple rows per position if multiple variants/categories)
cbpnet_pos_cat <- unique(
  chrompbnet_hg37_result[, .(chrpos, category)]
)

# All ChromBPNet positions (for diagnostics)
cbpnet_chrpos37 <- unique(cbpnet_pos_cat$chrpos)

autosomes <- as.character(1:22)

message("Total cbpnet chr:pos (hg19) with category info: ", nrow(cbpnet_pos_cat))
message("Example cbpnet chr:pos: ", paste(head(cbpnet_chrpos37), collapse = ", "))

## ----------------------------------------------------------------------
## 2) Helper: chr:pos from SuSiE IDs (e.g., '7:114700101:G:A' → '7:114700101')
## ----------------------------------------------------------------------

strip_to_chrpos <- function(x) sub(":[ACGT]+:[ACGT]+$", "", x)

## ----------------------------------------------------------------------
## 3) Region-level helper: CS-level membership + category
## ----------------------------------------------------------------------

membership_region <- function(region, gwas_id, cbpnet_pos_cat) {
  # If SuSiE didn’t converge or no credible sets, return empty
  if (isFALSE(region$converged) ||
      is.null(region$sets$cs) ||
      length(region$sets$cs) == 0) {
    return(data.table())
  }
  
  # SuSiE variant IDs like "7:114700101:G:A"
  susie_ids <- names(region$pip)
  if (is.null(susie_ids)) {
    stop("Region has no variant names in 'pip'.")
  }
  
  # Chromosome of this region
  region_chr <- unique(sub(":.*", "", susie_ids))
  if (length(region_chr) != 1L) {
    warning("Region has variants on multiple chromosomes? ", paste(region_chr, collapse = ", "))
  }
  region_chr <- region_chr[1]
  
  cs_list  <- region$sets$cs
  cs_names <- names(cs_list)
  if (is.null(cs_names)) {
    cs_names <- paste0("L", seq_along(cs_list))
  }
  
  # If region not on autosomes → CS exist but not eligible for CBPNet
  if (!(region_chr %in% autosomes)) {
    return(data.table(
      gwas_id    = gwas_id,
      chr        = region_chr,
      cs_id      = cs_names,
      eligible   = FALSE,
      has_cbpnet = FALSE,
      category   = NA_character_
    ))
  }
  
  # Region is autosomal → map SuSiE IDs to chr:pos
  susie_chrpos <- strip_to_chrpos(susie_ids)
  
  # Variant-level table
  dt_var <- data.table(
    idx    = seq_along(susie_ids),
    chrpos = susie_chrpos
  )
  
  # Attach ChromBPNet categories at positions where they exist
  # Result: one row per (variant index, category)
  dt_var_cat <- cbpnet_pos_cat[dt_var, on = "chrpos", nomatch = 0L, allow.cartesian = TRUE]
  # dt_var_cat has: chrpos, category, idx
  
  # Summarise at CS level
  cs_dt <- rbindlist(lapply(seq_along(cs_list), function(k) {
    idx_vec <- cs_list[[k]]
    
    cs_rows <- dt_var_cat[idx %in% idx_vec]
    has_cb  <- nrow(cs_rows) > 0
    
    cats <- unique(cs_rows[!is.na(category), category])
    
    cs_cat <- NA_character_
    if (has_cb) {
      if ("QTL_shared" %in% cats) {
        cs_cat <- "QTL_shared"
      } else if ("ChromBPNet_specific" %in% cats) {
        cs_cat <- "ChromBPNet_specific"
      } else {
        cs_cat <- NA_character_
      }
    }
    
    data.table(
      gwas_id    = gwas_id,
      chr        = region_chr,
      cs_id      = cs_names[k],
      eligible   = TRUE,
      has_cbpnet = has_cb,
      category   = cs_cat
    )
  }), use.names = TRUE, fill = TRUE)
  
  cs_dt
}

## ----------------------------------------------------------------------
## 4) Main: GWAS-level summary (per credible set, with categories)
## ----------------------------------------------------------------------

make_table <- function(gwas_id, cbpnet_pos_cat) {
  message("Processing GWAS: ", gwas_id)
  
  gwas_dir <- "/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/"
  
  # All region RDS files for this GWAS
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Restrict to the SuSiE 1.5Mb output
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  
  if (length(regions) == 0L) {
    warning("No SuSiE region files found for ", gwas_id)
    return(
      data.frame(
        gwas_id             = gwas_id,
        total_cs_all        = 0L,
        total_cs_eligible   = 0L,
        cs_with_cbpnet      = 0L,
        cs_cbpnet_specific  = 0L,
        cs_cbpnet_shared    = 0L,
        prop_cs_with_cbpnet = NA_real_
      )
    )
  }
  
  ## --- use the *same* CS counting as the original GWAS–QTL script ---
  count_cs_region <- function(file_name) {
    region <- readRDS(file_name)
    
    if(region$converged == FALSE){
      region_cs = 1
    } else {
      region_cs = length(region$sets$cs_index)
    }
    
    return(region_cs)
  }
  
  total_cs_all <- sum(sapply(regions, count_cs_region))
  ## -------------------------------------------------------------------
  ## -------------------------------------------------------------------
  
  ## --- GLOBAL POSITION-ONLY OVERLAP DIAGNOSTIC  ---
  all_susie_ids <- unique(unlist(lapply(regions, function(f) {
    region <- readRDS(f)
    names(region$pip)
  })))
  
  susie_chrpos_all <- unique(strip_to_chrpos(all_susie_ids))
  pos_overlap <- length(intersect(susie_chrpos_all, cbpnet_chrpos37))
  message("Position-only overlap (chr:pos) SuSiE ∩ cbpnet for ", gwas_id, ": ", pos_overlap)
  ## -----------------------------------------------------------
  
  # Per-CS table across regions 
  cs_tbl <- rbindlist(lapply(regions, function(f) {
    region <- readRDS(f)
    membership_region(region, gwas_id, cbpnet_pos_cat)
  }), use.names = TRUE, fill = TRUE)
  
  if (nrow(cs_tbl) == 0L) {
    return(
      data.frame(
        gwas_id             = gwas_id,
        total_cs_all        = total_cs_all,  # still report the SuSiE total
        total_cs_eligible   = 0L,
        cs_with_cbpnet      = 0L,
        cs_cbpnet_specific  = 0L,
        cs_cbpnet_shared    = 0L,
        prop_cs_with_cbpnet = NA_real_
      )
    )
  }
  
  total_cs_eligible <- sum(cs_tbl$eligible, na.rm = TRUE)
  cs_eligible <- cs_tbl[eligible == TRUE]
  
  n_cs_cbpnet   <- sum(cs_eligible$has_cbpnet)
  n_cs_specific <- sum(cs_eligible$has_cbpnet & cs_eligible$category == "ChromBPNet_specific", na.rm = TRUE)
  n_cs_shared   <- sum(cs_eligible$has_cbpnet & cs_eligible$category == "QTL_shared",           na.rm = TRUE)
  
  data.frame(
    gwas_id             = gwas_id,
    total_cs_all        = total_cs_all,      # now matches original total_cond_indep
    total_cs_eligible   = total_cs_eligible, # subset of those that are autosomal + analyzable
    cs_with_cbpnet      = n_cs_cbpnet,
    cs_cbpnet_specific  = n_cs_specific,
    cs_cbpnet_shared    = n_cs_shared,
    prop_cs_with_cbpnet = if (total_cs_eligible > 0) n_cs_cbpnet / total_cs_eligible else NA_real_
  )
}


## ----------------------------------------------------------------------
## 5) Run for one GWAS and then all GWAS
## ----------------------------------------------------------------------

# Quick test on a single GWAS
result_1 <- make_table(gwas_id, cbpnet_pos_cat)
print(result_1)

# Full table across all GWAS
all_results <- do.call(
  rbind,
  lapply(file_names_filt, make_table, cbpnet_pos_cat = cbpnet_pos_cat)
)
head(all_results)

fwrite(all_results,
       "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_immune_gwas_membership_all_results.csv")

all_results<-fread("~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_immune_gwas_membership_all_results.csv")
## ----------------------------------------------------------------------
## 6) Trait-level summary (with category split)
## ----------------------------------------------------------------------

trait_summary <- all_results %>%
  separate(gwas_id, into = c("year", "accession", "trait", "ethnicity"),
           sep = "_", remove = FALSE) %>%
  group_by(trait) %>%
  summarise(
    total_cs_all.mean        = mean(total_cs_all, na.rm = TRUE),
    total_cs_all.median      = median(total_cs_all, na.rm = TRUE),
    total_cs_all             = sum(total_cs_all, na.rm = TRUE),
    total_cs_eligible        = sum(total_cs_eligible, na.rm = TRUE),
    cs_with_cbpnet           = sum(cs_with_cbpnet, na.rm = TRUE),
    cs_cbpnet_specific       = sum(cs_cbpnet_specific, na.rm = TRUE),
    cs_cbpnet_shared         = sum(cs_cbpnet_shared, na.rm = TRUE),
    prop_cs_with_cbpnet_pooled = ifelse(
      total_cs_eligible > 0,
      cs_with_cbpnet / total_cs_eligible,
      NA_real_
    ),
    .groups = "drop"
  )

fwrite(trait_summary,
       "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_immune_gwas_membership.csv")

trait_summary <- fread("~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_immune_gwas_membership.csv")

# Keep only traits that exist in trait_labels
trait_summary <- trait_summary %>%
  filter(trait %in% names(trait_labels)) %>%
  mutate(
    trait_label = trait_labels[trait]
  )

# Ordering traits by total_cs_all.mean
trait_order <- trait_summary %>%
  distinct(trait, trait_label, total_cs_all.mean) %>%
  arrange(total_cs_all.mean) %>%
  pull(trait_label)

trait_summary <- trait_summary %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

## ----------------------------------------------------------------------
## 7) Long-format for stacked bars by category
## ----------------------------------------------------------------------

trait_long <- trait_summary %>%
  select(trait, trait_label, total_cs_eligible,
         cs_cbpnet_specific, cs_cbpnet_shared) %>%
  tidyr::pivot_longer(
    cols      = c(cs_cbpnet_specific, cs_cbpnet_shared),
    names_to  = "cbpnet_cat",
    values_to = "cs_n"
  ) %>%
  mutate(
    cbpnet_cat = dplyr::recode(cbpnet_cat,
                               cs_cbpnet_specific = "ChromBPNet_specific",
                               cs_cbpnet_shared   = "QTL_shared"
    ),
    prop = ifelse(total_cs_eligible > 0, cs_n / total_cs_eligible, NA_real_),
    trait_label = factor(trait_label, levels = trait_order)
  )

# For dots/line (average CS count per trait)
avg_total_df <- trait_summary %>%
  distinct(trait_label, total_cs_all.mean) %>%
  arrange(total_cs_all.mean) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

# Scaling between proportion and count axes
max_prop  <- max(trait_long$prop, na.rm = TRUE)
max_total <- max(avg_total_df$total_cs_all.mean, na.rm = TRUE)
scaling_factor <- max_prop / max_total

## ----------------------------------------------------------------------
## 8) Plot: stacked bars by category + dots/line for total CS
## ----------------------------------------------------------------------

scaling_factor <- 1 / 500

plot2 <- ggplot(trait_long, aes(x = trait_label, y = prop, fill = cbpnet_cat)) +
  geom_col(color = "gray40", size = 0.1) +
  geom_text(
    data = trait_summary,
    aes(x = trait_label,
        y = prop_cs_with_cbpnet_pooled / 2,
        label = cs_with_cbpnet),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  geom_point(
    data = avg_total_df,
    aes(x = trait_label, y = total_cs_all.mean * scaling_factor),
    inherit.aes = FALSE,
    shape = 21,
    fill = "gray30",
    size = 3,
    stroke = 0.3,
    alpha = 0.5
  ) +
  geom_line(
    data = avg_total_df,
    aes(x = trait_label, y = total_cs_all.mean * scaling_factor, group = 1),
    inherit.aes = FALSE,
    color = "gray50",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  scale_fill_manual(
    values = c(
      "ChromBPNet_specific" = "#0071ba",
      "QTL_shared"          = "#ff9f1c"
    ),
    name = "ChromBPNet category"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    name = "Proportion of GWAS credible sets with CBPNet variants",
    sec.axis = sec_axis(
      ~ . / scaling_factor,
      # choose either of these:
      # breaks = seq(0, 500, 100),
      breaks = seq(100, 500, 100),   # if you only want 100–500 shown
      name   = "Average number of credible sets per trait"
    )
  ) +
  coord_flip() +
  labs(
    title = "ChromBPNet variant membership in immune GWAS credible sets",
    x = "GWAS Trait"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x        = element_text(hjust = 1, size = 16),
    axis.text.y        = element_text(size = 16),
    axis.title         = element_text(size = 16, face = "bold"),
    legend.position    = "right",
    legend.justification = "center",
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 12),
    axis.text.x.bottom = element_text(size = 14, color = "black"),
    axis.text.x.top    = element_text(size = 12, color = "gray50"),
    axis.title.x.top   = element_text(size = 14, color = "gray50")
  )

print(plot2)


ggsave(
  "~/cd4_qtl_paper_figures/figure_6/plotting/plots/cbpnet_gwas_coloc_proportions_by_category.pdf",
  plot2, height = 15, width = 8, dpi = 300, units = "in"
)




