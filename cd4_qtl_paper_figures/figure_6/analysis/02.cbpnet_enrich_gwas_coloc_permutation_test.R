.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
options(bitmapType = "cairo")
library(ggplot2)
library(data.table)
library(tidyverse)
library(forcats)
library(dplyr)
library(GenomicRanges)
library(lmtest)
library(sandwich)
library(future.apply)
library(future)
library(tibble)
library(tidyr)
library(ggplot2)
library(grid) 
n <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
plan(multisession, workers = max(1L, n))

source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")



count_cs = function(file_name){
  # We need the total conditionally independent GWAS variants
  # Load region rds file
  region = readRDS(file_name)

  if(region$converged == F){
    region_cs = 1
  } else {
    region_cs = length(region$sets$cs_index)
  }

  return(region_cs)
}


make_table <- function(gwas_id) {
  message("Processing: ", gwas_id)
  gwas_dir <- '/gcgnl/finemapping_autoimmune/'
  # Get region RDS files
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  gwas_cs <- sapply(regions, count_cs)

  # Check if coloc result files exist
  eqtl_path  <- paste0(data.dir, "coloc_results/eqtl_gwas_coloc/", gwas_id, "_coloc_results.csv")
  caqtl_path <- paste0(data.dir, "coloc_results/CD4T_chromatin_", gwas_id, "_coloc_results.csv")
  if (!file.exists(eqtl_path) || !file.exists(caqtl_path)) return(NULL)

  # Read coloc results
  eqtl_coloc_result <- fread(eqtl_path) %>%
    mutate(region_cs = paste(region, idx1, sep = "_")) %>%
    filter(pval < 1e-5, pval_nominal < 1e-3) %>%
    group_by(eQTL_variant_GRC37, gene) %>%                 # removing duplicates in unresolved regions
    slice_max(order_by = PP.H4.abf, n = 1, with_ties = FALSE) %>%  # keep highest PP.H4.abf
    ungroup() %>%
    mutate(gwas_id = gwas_id)

  caqtl_coloc_result <- fread(caqtl_path) %>%
    mutate(region_cs = paste(region, idx1, sep = "_")) %>%
    filter(pval < 1e-5, pval_nominal < 1e-3) %>%
    group_by(eQTL_variant_GRC37, peak) %>%                 # removing duplicates in unresolved regions
    slice_max(order_by = PP.H4.abf, n = 1, with_ties = FALSE) %>%  # keep highest PP.H4.abf
    ungroup() %>%
    mutate(gwas_id = gwas_id)

  shared_coloc_result <- caqtl_coloc_result %>%
    filter(region_cs %in% eqtl_coloc_result$region_cs)  %>%
    group_by(eQTL_variant_GRC37, peak) %>%                 # removing duplicates in unresolved regions
    slice_max(order_by = PP.H4.abf, n = 1, with_ties = FALSE) %>%  # keep highest PP.H4.abf
    ungroup()%>%
    mutate(gwas_id = gwas_id)

  #get unique sets of regions
  only_eqtl <- setdiff(unique(eqtl_coloc_result$region_cs), caqtl_coloc_result$region_cs)
  only_caqtl <- setdiff(unique(caqtl_coloc_result$region_cs), eqtl_coloc_result$region_cs)
  both <- unique(shared_coloc_result$region_cs)

  #filter back
  only_eqtl_df <- eqtl_coloc_result[eqtl_coloc_result$region_cs %in% only_eqtl, ]
  only_caqtl_df <- caqtl_coloc_result[caqtl_coloc_result$region_cs %in% only_caqtl, ]
  both_df <- shared_coloc_result[shared_coloc_result$region_cs %in% both, ]

  intersect(only_eqtl, both)         # should be length 0
  intersect(only_caqtl, both)        # should be length 0
  intersect(only_eqtl, only_caqtl)   # should be length 0

  only_eqtl_df <- only_eqtl_df %>%
    select(GWAS_variant_GRC37, region_cs, variant_id_GRC38, region_GRC38, idx2, gene, gwas_id)

  only_caqtl_df <- only_caqtl_df %>%
    select(GWAS_variant_GRC37, region_cs, variant_id_GRC38, region_GRC38, idx2, peak, gwas_id)

  both_df <- both_df %>%
    select(GWAS_variant_GRC37, region_cs, variant_id_GRC38, region_GRC38, idx2, peak, gwas_id)


  return(list(eqtl_only_ls=only_eqtl_df, only_caqtl_ls=only_caqtl_df, both_ls=both_df ))
}


data.dir="/gpfs/commons/groups/lappalainen_lab/sghatan/marlis_pj/coloc/"
# folder names
file_names = unlist(fread(paste0(data.dir, "preprocessed_folder_names.txt"), header = F))
#gwas_id = head(file_names,20) #for testing purposes



variant_list <- lapply(file_names, make_table)

results_filtered <- variant_list[!sapply(variant_list, function(x) {
  is.null(x) || (is.atomic(x) && any(is.na(x))) || (is.data.frame(x) && ncol(x) == 0)
})]

# Combine and save summary table
eqtl_only_df <- bind_rows(lapply(results_filtered, `[[`, "eqtl_only_ls")) %>%
  mutate(finemapped_cs_eqtl = paste0(gene, "_",  region_GRC38, "_", idx2))
fwrite(eqtl_only_df, "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_eqtl_only_df.csv")


# Combine and save shared coloc results
only_caqtl_df <- bind_rows(lapply(results_filtered, `[[`, "only_caqtl_ls")) %>%
  mutate(finemapped_cs_caqtl = paste0(peak, "_",  region_GRC38, "_", idx2))
fwrite(only_caqtl_df, "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_caqtl_only_df.csv")


both <- bind_rows(lapply(results_filtered, `[[`, "both_ls")) %>%
  mutate(finemapped_cs_caqtl = paste0(peak, "_",  region_GRC38, "_", idx2))
fwrite(both, "~/cd4_qtl_paper_figures/figure_6/data/CD4T_coloc_summary_table_both_caqtl_eqtl_df.csv")


##map lead variant to finemapping results
finemap <- fread("/gchm/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv",
                               data.table = FALSE, stringsAsFactors = FALSE)


##please note that regions for the GWAS-colocalization are different from the regions of the finemapping, so keeping track of the right columns is important

## get the variants for each credible set.
only_eqtl2_finemap <- eqtl_only_df %>%
  left_join(
    finemap %>% select(peak, variant_id.y, variant_id.x, eQTL_variant,caQTL_variant, chromosome,  variant_pos.y, variant_pos.x,  finemapped_cs_eqtl, finemapped_cs_coloc ),   #variant_id.y are all the finemapped variants
    by = "finemapped_cs_eqtl"
  ) %>% mutate(coloc_type="eQTL_only")

colnames(only_eqtl2_finemap) <-c("GWAS_variant_GRC37" , "region_cs"  ,    "variant_id_GRC38" ,  "region_GRC38" ,  "idx2"  ,   "feature_1",    "gwas_id"   ,    "finemapped_cs_qtl" , "feature_2"     ,   "variant_id.y" ,
  "variant_id.x"    ,    "eQTL_variant"  ,   "caQTL_variant" ,   "chromosome"    ,    "variant_pos.y" ,  "variant_pos.x" ,   "finemapped_cs_coloc", "coloc_type")

print(paste0("Number of unique eQTL only variants before merging finemapping data: ", length(unique(eqtl_only_df$finemapped_cs_eqtl))))
print(paste0("Number of unique eQTL only after merging finemapping data: ", length(unique(only_eqtl2_finemap$finemapped_cs_qtl))))


only_caqtl2_finemap <- only_caqtl_df %>%
  left_join(
    finemap %>% select(gene, variant_id.y, variant_id.x, eQTL_variant, caQTL_variant,  chromosome,   variant_pos.y, variant_pos.x,  finemapped_cs_caqtl, finemapped_cs_coloc),   #variant_id.y are all the finemapped variants
    by = "finemapped_cs_caqtl"
  ) %>% mutate(coloc_type="caQTL_only")

colnames(only_caqtl2_finemap)<-c("GWAS_variant_GRC37" , "region_cs"  ,    "variant_id_GRC38" ,  "region_GRC38" ,  "idx2"  ,   "feature_1",    "gwas_id"   ,    "finemapped_cs_qtl" , "feature_2"     ,   "variant_id.y" ,
                                 "variant_id.x"    ,    "eQTL_variant"  ,   "caQTL_variant" ,   "chromosome"    ,    "variant_pos.y" ,  "variant_pos.x" ,   "finemapped_cs_coloc", "coloc_type")

print(paste0("Number of unique caQTL only variants before merging finemapping data: ", length(unique(only_caqtl_df$finemapped_cs_caqtl))))
print(paste0("Number of unique caQTL only after merging finemapping data: ",length(unique(only_caqtl2_finemap$finemapped_cs_qtl))) )


both2_finemap <- both  %>%
  left_join(
    finemap %>% select(gene, variant_id.y, variant_id.x, eQTL_variant, caQTL_variant,  chromosome,  variant_pos.y, variant_pos.x, finemapped_cs_caqtl, finemapped_cs_coloc),   #variant_id.y are all the finemapped variants
    by = "finemapped_cs_caqtl") %>% mutate(coloc_type="both")

colnames(both2_finemap) <-c("GWAS_variant_GRC37" , "region_cs"  ,    "variant_id_GRC38" ,  "region_GRC38" ,  "idx2"  ,   "feature_1",    "gwas_id"   ,    "finemapped_cs_qtl" , "feature_2"     ,   "variant_id.y" ,
                            "variant_id.x"    ,    "eQTL_variant"  ,   "caQTL_variant" ,   "chromosome"    ,    "variant_pos.y" ,  "variant_pos.x" ,   "finemapped_cs_coloc", "coloc_type")

print(paste0("Number of unique variants colocalized for both eQTL and caQTL only before merging finemapping data: ", length(unique(both$finemapped_cs_caqtl))))
print(paste0("Number of unique variants colocalized for both eQTL and caQTL after merging finemapping data: ", length(unique(both2_finemap$finemapped_cs_qtl))))


## my goal was to match the lead variants here to a credible set, so i can get the variants in the credible set and then see if any of these variants are chrombpnet variants,
## then i would get summarize the same based on unique gwas variants w/ or w/o chrombpnet variants

print("Reading ChromBPNet Variants")
##Read chromBPnet variants
cbpnet <- fread('/gchm/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores_significan.tsv',
                data.table = FALSE, stringsAsFactors = FALSE)

cbpnet_hits  <- unique((cbpnet[cbpnet$sig_vars_IPS_p0.05==TRUE,])$variant_id)
cbpnet_callable <- unique(cbpnet$variant_id)

print(paste0("There are ", length(cbpnet_callable), " ChromBPNet callable variants"))
print(paste0("There are ", length(cbpnet_hits), " ChromBPNet significantly scored variants"))

#get the variants for each set
find_cbpnet_hits <- function(df, cbpnet_hits) {
  df %>%
    select(variant_id.y, variant_id.x, eQTL_variant, caQTL_variant) %>%
    pivot_longer(cols = everything(), values_to = "variant_id") %>%
    filter(!is.na(variant_id)) %>%
    distinct(variant_id) %>%
    mutate(
      variant_id_normalized = gsub(
        "([0-9]+):([0-9]+)\\[b38\\]([A-Z]+),([A-Z]+)",
        "chr\\1_\\2_\\3_\\4",
        variant_id
      ),
      chrombpnet_hit = variant_id_normalized %in% cbpnet_hits,
      chrombpnet_callable = variant_id_normalized %in% cbpnet_callable
    )
}

only_eqtl_vars   <- find_cbpnet_hits(only_eqtl2_finemap, cbpnet_hits)
only_caqtl_vars  <- find_cbpnet_hits(only_caqtl2_finemap, cbpnet_hits)
both2_finemap_vars <- find_cbpnet_hits(both2_finemap,        cbpnet_hits)


print("Annotating colocalized variants for ChromBPNet hits (true/false)")
# 1) ID normalizer (once)
normalize_id <- function(x) {
  x <- as.character(x)
  sub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]+),([A-Z]+)", "chr\\1_\\2_\\3_\\4", x, perl = TRUE)
}

# 2) Quick name sanitizer (fixes the dplyr error)
sanitize_names <- function(df) {
  nm <- names(df)
  bad <- is.na(nm) | nm == ""
  if (any(bad)) names(df)[bad] <- paste0(".col", seq_len(sum(bad)))
  df
}

# 3) Vectorized annotator that takes hits_df (your usage)
annotate_original_df_fast_hitsdf <- function(original_df, hits_df,
                                             id_cols = c("variant_id.y", "variant_id.x",
                                                         "eQTL_variant", "caQTL_variant")) {
  original_df <- sanitize_names(original_df)
  stopifnot(all(id_cols %in% names(original_df)))

  hit_ids <- hits_df %>%
    filter(chrombpnet_hit) %>%
    pull(variant_id_normalized) %>%
    unique()

  callable_ids <- hits_df %>%
    filter(chrombpnet_callable) %>%
    pull(variant_id_normalized) %>%
    unique()

  norm_cols <- paste0(id_cols, "__norm")

  original_df %>%
    mutate(
      across(all_of(id_cols), ~ normalize_id(.x), .names = "{.col}__norm")
    ) %>%
    mutate(
      chrombpnet_hit_any      = if_any(all_of(norm_cols), ~ replace_na(.x, "") %in% hit_ids),
      chrombpnet_callable_any = if_any(all_of(norm_cols), ~ replace_na(.x, "") %in% callable_ids)
    ) %>%
    select(-all_of(norm_cols))
}


only_eqtl2_finemap_annotated  <- annotate_original_df_fast_hitsdf(only_eqtl2_finemap,  only_eqtl_vars)
only_caqtl2_finemap_annotated <- annotate_original_df_fast_hitsdf(only_caqtl2_finemap, only_caqtl_vars)
both2_finemap_annotated       <- annotate_original_df_fast_hitsdf(both2_finemap,       both2_finemap_vars)



variants_merged <- bind_rows(only_eqtl2_finemap_annotated,
                             only_caqtl2_finemap_annotated,
                             both2_finemap_annotated)


# Ensure key CBPNet columns are logical
variants_merged <- variants_merged %>%
  mutate(
    chrombpnet_callable_any = as.logical(chrombpnet_callable_any),
    chrombpnet_hit_any      = as.logical(chrombpnet_hit_any)
  )


# Exclude extended MHC (optional but recommended)
not_mhc <- !(variants_merged$chromosome %in% c("chr6") &
               variants_merged$`variant_pos.y` >= 25e6 &
               variants_merged$`variant_pos.x` >= 25e6 &
               variants_merged$`variant_pos.y` <= 34e6 &
               variants_merged$`variant_pos.x` <= 34e6)

vm <- variants_merged[not_mhc, ]

# CS-level aggregation using *your* columns
# ---- 2) Build CS-level summaries for the three categories ----
cs_coloc <- vm %>%
  group_by(finemapped_cs_coloc, coloc_type, region_GRC38, chromosome) %>%
  summarise(
    m_callable = sum(chrombpnet_callable_any, na.rm = TRUE),
    h_hits     = sum(chrombpnet_hit_any,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(frac_hits = ifelse(m_callable > 0, h_hits / m_callable, NA_real_))
# Quick sanity: # hits should be ≤ # callable
stopifnot(all(cs_coloc$h_hits <= cs_coloc$m_callable, na.rm = TRUE))

###### SETTING UP THE BACKGROUND
##non-colocalized credible sets
non_colocalized <- finemap %>% filter(!(finemapped_cs_coloc %in% variants_merged$finemapped_cs_coloc)) %>% select( "region.x" ,  "eQTL_variant"  ,   "caQTL_variant" ,  "variant_id.y", "variant_id.x",  "chromosome" ,    "variant_pos.y" ,  "variant_pos.x" ,   "finemapped_cs_coloc")
non_colocalized_vars <- find_cbpnet_hits(non_colocalized,        cbpnet_hits)
non_colocalized_annotated <- annotate_original_df_fast_hitsdf(non_colocalized,  non_colocalized_vars)

baseline_vm <- non_colocalized_annotated %>%
  mutate(
    chrombpnet_callable_any = as.logical(chrombpnet_callable_any),
    chrombpnet_hit_any      = as.logical(chrombpnet_hit_any)
  )

# same xMHC exclusion
not_mhc_b <- !(baseline_vm$chromosome %in% c("chr6") &
                 baseline_vm$`variant_pos.y` >= 25e6 &
                 baseline_vm$`variant_pos.x` >= 25e6 &
                 baseline_vm$`variant_pos.y` <= 34e6 &
                 baseline_vm$`variant_pos.x` <= 34e6)
baseline_vm <- baseline_vm[not_mhc_b, ]

baseline_cs <- baseline_vm %>%
  group_by(finemapped_cs_coloc, region.x , chromosome) %>%
  summarise(
    coloc_type ="non_coloc",
    region_GRC38 =region.x,
    m_callable = sum(chrombpnet_callable_any, na.rm = TRUE),
    h_hits     = sum(chrombpnet_hit_any,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(frac_hits = ifelse(m_callable > 0, h_hits / m_callable, NA_real_)) %>% select("finemapped_cs_coloc","coloc_type", "region_GRC38", "chromosome", "m_callable", "h_hits", "frac_hits")

#same columns names
colnames(baseline_cs)==colnames(cs_coloc)

write_csv(cs_coloc, "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_var_gwas_enrichment_glm_coloc_test.csv")
write_csv(baseline_cs, "~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_var_gwas_enrichment_glm_coloc_background.csv")

cs_coloc<-fread("~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_var_gwas_enrichment_glm_coloc_test.csv")
baseline_cs<-fread("~/cd4_qtl_paper_figures/figure_6/data/chrombpnet_var_gwas_enrichment_glm_coloc_background.csv")

##Enrichment
print("Starting enrichment for chromBPNet variants in GWAS colocalized QTL loci")


# -------- helpers --------
call_bin_fun <- function(x) cut(x, c(1,2,5,10,20,50,Inf), include.lowest=TRUE, right=FALSE)

# -------- test one category vs baseline (GLM + binomial simulation null) --------
test_category_fast <- function(cs_coloc, baseline_cs, category_label,
                               nperm = 10000L, seed = 1L,
                               save_rds = NULL) {
  # ---- Stage A: filter & prep ----
  A <- cs_coloc %>%
    dplyr::filter(coloc_type == category_label, m_callable >= 1) %>%
    dplyr::mutate(chromosome = as.character(chromosome),
                  frac_hits  = h_hits / m_callable)
  
  B <- baseline_cs %>%
    dplyr::filter(m_callable >= 1) %>%
    dplyr::mutate(chromosome = as.character(chromosome),
                  frac_hits  = h_hits / m_callable)
  
  stopifnot(nrow(A) > 0, nrow(B) > 0)
  
  # ---- Stage B: GLM (same as before) ----
  AB <- dplyr::bind_rows(dplyr::mutate(A, group="cat"),
                         dplyr::mutate(B, group="base"))
  
  m_comp <- stats::glm(
    cbind(h_hits, pmax(0, m_callable - h_hits)) ~ I(group=="cat") +
      factor(chromosome) + log(m_callable),
    family = binomial(), data = AB
  )
  naive <- summary(m_comp)$coefficients['I(group == "cat")TRUE', ]
  
  # ---- Stage C: permutation null (same logic) ----
  p_base_chr <- B %>%
    dplyr::group_by(chromosome) %>%
    dplyr::summarise(H = sum(h_hits), M = sum(m_callable), .groups="drop") %>%
    dplyr::mutate(p = ifelse(M > 0, H / M, 0)) %>%
    { stats::setNames(.$p, .$chromosome) }
  
  p_global <- sum(B$h_hits) / sum(B$m_callable)
  p_vec <- unname(ifelse(A$chromosome %in% names(p_base_chr),
                         p_base_chr[A$chromosome], p_global))
  
  set.seed(seed)
  m_j <- A$m_callable
  K   <- as.integer(nperm)
  
  # Vectorized simulation: generate K means of h*/m
  sim_means <- replicate(K, {
    h_sim <- stats::rbinom(length(m_j), size = m_j, prob = p_vec)
    mean(h_sim / m_j)
  })
  
  obs_mean <- mean(A$frac_hits)
  p_emp_one_sided <- (sum(sim_means >= obs_mean) + 1) / (K + 1)
  perm_null_mean  <- mean(sim_means)
  delta_frac      <- obs_mean - perm_null_mean
  
  # ---- Return: summary + plot-ready payload ----
  out <- tibble::tibble(
    category       = category_label,
    n_cs           = nrow(A),
    total_callable = sum(A$m_callable),
    total_hits     = sum(A$h_hits),
    mean_frac      = obs_mean,
    glm_logOR      = unname(naive["Estimate"]),
    glm_OR         = exp(unname(naive["Estimate"])),
    glm_CI_low     = exp(unname(naive["Estimate"] - 1.96 * naive["Std. Error"])),
    glm_CI_high    = exp(unname(naive["Estimate"] + 1.96 * naive["Std. Error"])),
    glm_p_two      = unname(naive["Pr(>|z|)"]),
    perm_null_mean = perm_null_mean,
    perm_delta     = delta_frac,
    perm_p_one     = p_emp_one_sided,
    # New: plot payload
    observed_stat  = obs_mean,
    null_stats     = list(sim_means),
    nperm          = K,
    seed           = seed
  )
  
  if (!is.null(save_rds)) {
    saveRDS(out, save_rds)
  }
  out
}


# ----------  Run all three categories and save results  ----------
# Run your three tests
res_eqtl  <- test_category_fast(cs_coloc, baseline_cs,category_label= "eQTL_only",  nperm = 1e5, seed = 11)
res_caqtl <- test_category_fast(cs_coloc, baseline_cs,category_label= "caQTL_only", nperm = 1e5, seed = 22)
res_both  <- test_category_fast(cs_coloc, baseline_cs, category_label="both",       nperm = 1e5, seed = 33)

# Combine results
all_res <- dplyr::bind_rows(res_both, res_eqtl, res_caqtl)
# Save to disk (RDS)
saveRDS(all_res,
        file = "~/cd4_qtl_paper_figures/figure_6/data/cbpnet_enrichment_results_by_category.rds")

all_res<-readRDS("~/cd4_qtl_paper_figures/figure_6/data/cbpnet_enrichment_results_by_category.rds")



colors <- c(
  eQTL_only  = "#c7197c",
  caQTL_only = "#9ccb86",
  both       = "#f39c12"
)

# ----- Build plotting DF (ASCII-safe labels) -----
df_or <- all_res %>%
  mutate(
    category = factor(category, levels = c("both", "eQTL_only", "caQTL_only")),
    perm_p_lab = ifelse(
      perm_p_one <= 1/nperm,
      paste0("perm p < ", format(1/nperm, scientific = TRUE)),
      paste0("perm p = ", signif(perm_p_one, 3))
    ),
    lab = sprintf(
      "OR = %.2f (%.2f-%.2f); %s; obs = %.3f",
      glm_OR, glm_CI_low, glm_CI_high, perm_p_lab, observed_stat
    )
  )

# constant x position for all labels (just right of the widest CI)
x_text <- max(df_or$glm_CI_high, na.rm = TRUE) * 1.08

# ----- Plot -----
p_or <- ggplot(df_or, aes(y = category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  
  geom_errorbarh(
    aes(xmin = glm_CI_low, xmax = glm_CI_high, color = category),
    height = 0.18,
    linewidth = 0.5
  ) +
  geom_point(aes(x = glm_OR, color = category), size = 5) +
  
  # label as constant x (NOT inside aes)
  geom_text(
    aes(label = lab),
    x = x_text,
    hjust = 1,
    vjust = -2,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = colors) +
  coord_cartesian(
    xlim = c(min(df_or$glm_CI_low, na.rm = TRUE) * 0.95, x_text * 1.02),
    clip = "off"
  ) +
  labs(
    x = "Odds ratio (GLM) with 95% CI",
    y = NULL,
    title = "ChromBPNet variant enrichment across colocalization classes"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(5.5, 80, 5.5, 5.5)  # extra right margin for labels
  )

p_or


ggsave(
  "~/cd4_qtl_paper_figures/figure_4/plotting/plots/chrombpnet_qtl_glm_OR_forest.pdf",
  p_or, width = 7, height = 3.5
)




# df_or <- all_res %>%
#   mutate(
#     category = factor(category, levels = c("both", "eQTL_only", "caQTL_only")),
#     lab = sprintf("OR %.2f (%.2f–%.2f)", glm_OR, glm_CI_low, glm_CI_high)
#   )
# 
# colors <- c(eQTL_only="#c7197c", caQTL_only="#9ccb86", both="#f39c12")
# x_text <- max(df_or$glm_CI_high) * 1.25
# 
# p_lolli <- ggplot(df_or, aes(y = category)) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
#   geom_segment(aes(x = 1, xend = glm_OR, yend = category),
#                linewidth = 1.2, color = "black", lineend = "round") +
#   geom_errorbarh(aes(xmin = glm_CI_low, xmax = glm_CI_high, color = category),
#                  height = 0.18, linewidth = 0.9) +
#   geom_point(aes(x = glm_OR, color = category), size = 6) +
#   geom_text(aes(x = x_text, label = lab), hjust = 0, size = 4, color = "black") +
#   scale_color_manual(values = colors) +
#   scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25","0.5","1","2","4")) +
#   coord_cartesian(xlim = c(min(df_or$glm_CI_low) * 0.9, x_text * 1.02)) +
#   labs(x = "Odds ratio (log scale)", y = NULL) +
#   theme_classic(base_size = 14) +
#   theme(legend.position = "none")
# 
# p_lolli
# 
# 
# ggsave("~/cd4_qtl_paper_figures/figure_4/plotting/plots/chrombpnet_qtl_perm_enrichment_lollipop.pdf",  lp_enrich, width = 5, height = 4)
# 
