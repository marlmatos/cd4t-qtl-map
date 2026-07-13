#!/usr/bin/env Rscript

# ======================================================================
# STANDALONE ANALYSIS:
# ChromBPNet prioritization enrichment in molQTL-colocalized GWAS loci
#
# PRIMARY QUESTION
# Among ChromBPNet-callable variants contained in GWAS credible sets,
# are ChromBPNet-prioritized variants enriched or depleted according to
# molecular-QTL colocalization category?
#
# Categories:
#   both        = eQTL and caQTL colocalization
#   just_caqtl  = caQTL colocalization only
#   just_eqtl   = eQTL colocalization only
#   no_coloc    = no eQTL or caQTL colocalization
#
# IMPORTANT DESIGN CHOICE
# The numerator and denominator are defined using GWAS credible-set
# variants only for every category. 
#
# Outputs:
#   1. GWAS_CS_only_callable_hit_counts_by_locus.csv
#   2. GWAS_CS_only_category_summary.csv
#   3. GWAS_CS_only_correct_enrichment_results.csv
#   4. GWAS_CS_only_quasibinomial_model.rds
#   5. GWAS_CS_only_permutation_results_with_nulls.rds
#   6. GWAS_CS_only_global_colocalization_factor_test.txt
#   7. GWAS_CS_only_correct_enrichment_forest_plot.pdf/png
#   8. Multiple QC tables
# ======================================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
options(bitmapType = "cairo")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(scales)
})

# ======================================================================
# 0) CONFIGURATION: EDIT THESE PATHS
# ======================================================================
source("~/cd4t-qtl-map/analysis/paths_config.R")

data_dir<-script13$data_dir
gwas_fm_dir<-script13$gwas_fm_dir
chrombpnet_score_file<-script13$chrombpnet_score_file
chrombpnet_mapping_file<-script13$chrombpnet_mapping_file
eqtl_cs_dir<-script13$eqtl_cs_dir
caqtl_cs_dir<-script13$caqtl_cs_dir


output_dir <- file.path(
  data_dir,
  "ChromBPNet_GWAS_CS_only_enrichment"
)

# Name of the Boolean or Boolean-convertible column defining
# ChromBPNet-prioritized variants.
CHROMBPNET_HIT_COLUMN <- "sig_vars_IPS_p0.05"

GWAS_P_THRESH <- 1e-5
QTL_P_THRESH <- 1e-3
PP_H4_THRESH <- 0.5

N_PERM <- 100000L

PERM_SEEDS <- c(
  both = 11L,
  just_caqtl = 22L,
  just_eqtl = 33L
)

CALLABLE_BREAKS <- c(1, 2, 5, 10, 20, 50, Inf)

EXCLUDE_EXTENDED_MHC <- TRUE
MHC_CHR <- "6"
MHC_START <- 25e6
MHC_END <- 34e6

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

configured_paths <- c(
  data_dir,
  gwas_fm_dir,
  chrombpnet_score_file,
  chrombpnet_mapping_file,
  eqtl_cs_dir,
  caqtl_cs_dir
)

if (any(grepl("REPLACE_WITH", configured_paths, fixed = TRUE))) {
  stop(
    "Edit every REPLACE_WITH entry in the CONFIGURATION section before running."
  )
}

required_files <- c(
  chrombpnet_score_file,
  chrombpnet_mapping_file,
  file.path(data_dir, "preprocessed_folder_names.txt")
)

required_dirs <- c(
  data_dir,
  gwas_fm_dir,
  eqtl_cs_dir,
  caqtl_cs_dir
)

missing_files <- required_files[!file.exists(required_files)]
missing_dirs <- required_dirs[!dir.exists(required_dirs)]

if (length(missing_files) > 0L) {
  stop(
    "Missing input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

if (length(missing_dirs) > 0L) {
  stop(
    "Missing input directory/directories: ",
    paste(missing_dirs, collapse = ", ")
  )
}

# ======================================================================
# 1) CONSTANT KEYS AND GENERAL HELPERS
# ======================================================================

KEY_STUDY <- c(
  "gwas_id",
  "trait",
  "region",
  "cs_index",
  "region_cs"
)

KEY_LOCUS <- c(
  "trait",
  "region",
  "cs_index",
  "region_cs"
)

KEY_QTL <- c(
  "qtl_type",
  "feature",
  "qtl_region",
  "qtl_cs_index"
)

extract_trait <- function(gwas_id) {
  x <- sub("_preprocessed$", "", gwas_id)
  pieces <- strsplit(x, "_", fixed = TRUE)[[1]]
  
  if (length(pieces) < 3L) {
    stop("Cannot parse trait from GWAS ID: ", gwas_id)
  }
  
  pieces[3]
}

make_region_cs <- function(region, cs_index) {
  paste(region, cs_index, sep = "_")
}

collapse_unique_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  x <- sort(unique(x))
  
  if (length(x) == 0L) {
    NA_character_
  } else {
    paste(x, collapse = ";")
  }
}

normalize_hg19_key <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x)
  
  x <- sub(
    "^([^_]+)_([0-9]+)_([^_]+)_([^_]+)$",
    "\\1:\\2:\\3:\\4",
    x
  )
  
  x
}

normalize_hg38_key <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x)
  
  x <- sub(
    "^([^:]+):([0-9]+)\\[b38\\]([^,]+),(.+)$",
    "\\1:\\2:\\3:\\4",
    x
  )
  
  x <- sub(
    "^([^_]+)_([0-9]+)_([^_]+)_([^_]+)$",
    "\\1:\\2:\\3:\\4",
    x
  )
  
  x
}

parse_variant_chr <- function(x) {
  sub(":.*$", "", normalize_hg19_key(x))
}

parse_variant_pos <- function(x) {
  key <- normalize_hg19_key(x)
  
  valid <- grepl("^[^:]+:[0-9]+:", key)
  
  out <- rep(NA_real_, length(key))
  out[valid] <- suppressWarnings(
    as.numeric(
      sub(
        "^[^:]+:([0-9]+):.*$",
        "\\1",
        key[valid]
      )
    )
  )
  
  out
}

make_callable_bin <- function(x) {
  cut(
    x,
    breaks = CALLABLE_BREAKS,
    include.lowest = TRUE,
    right = FALSE
  )
}

safe_fraction <- function(h, m) {
  denom <- sum(m, na.rm = TRUE)
  
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  
  sum(h, na.rm = TRUE) / denom
}

load_cs_files <- function(dir, pattern) {
  files <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0L) {
    stop(
      "No CS files found in: ",
      dir,
      " using pattern: ",
      pattern
    )
  }
  
  message("Loading ", length(files), " CS files from: ", dir)
  
  out <- data.table::rbindlist(
    lapply(files, data.table::fread),
    fill = TRUE
  )
  
  if ("cs" %in% names(out) && !"cs_index" %in% names(out)) {
    data.table::setnames(out, "cs", "cs_index")
  }
  
  out
}

# ======================================================================
# 2) LOAD AND CLEAN CHROMBPNET RESULTS
# ======================================================================

cbpnet_scores <- fread(
  chrombpnet_score_file,
  data.table = FALSE
)

cbpnet_mapping <- fread(
  chrombpnet_mapping_file,
  data.table = FALSE
)

required_score_cols <- c(
  "variant_id",
  CHROMBPNET_HIT_COLUMN
)

required_mapping_cols <- c(
  "variant_id",
  "hg19_chr",
  "hg19_pos",
  "alt",
  "ref"
)

missing_score_cols <- setdiff(
  required_score_cols,
  names(cbpnet_scores)
)

missing_mapping_cols <- setdiff(
  required_mapping_cols,
  names(cbpnet_mapping)
)

if (length(missing_score_cols) > 0L) {
  stop(
    "ChromBPNet score file is missing: ",
    paste(missing_score_cols, collapse = ", ")
  )
}

if (length(missing_mapping_cols) > 0L) {
  stop(
    "ChromBPNet mapping file is missing: ",
    paste(missing_mapping_cols, collapse = ", ")
  )
}

cbpnet_scores_clean <- cbpnet_scores %>%
  transmute(
    variant_id = as.character(variant_id),
    chrombpnet_hit = as.logical(
      .data[[CHROMBPNET_HIT_COLUMN]]
    )
  ) %>%
  mutate(
    chrombpnet_hit = replace_na(
      chrombpnet_hit,
      FALSE
    )
  ) %>%
  group_by(variant_id) %>%
  summarise(
    chrombpnet_hit = any(chrombpnet_hit),
    .groups = "drop"
  )

cbpnet_mapping_clean <- cbpnet_mapping %>%
  transmute(
    variant_id = as.character(variant_id),
    hg19_chr = gsub(
      "^chr",
      "",
      as.character(hg19_chr)
    ),
    hg19_var_lookup = paste0(
      hg19_chr,
      ":",
      hg19_pos,
      ":",
      alt,
      ":",
      ref
    ),
    variant_id_b38 = sub(
      "^chr([^_]+)_([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$",
      "\\1:\\2[b38]\\3,\\4",
      variant_id
    )
  ) %>%
  filter(
    !is.na(variant_id),
    variant_id != "",
    !is.na(hg19_var_lookup),
    hg19_var_lookup != "",
    !is.na(variant_id_b38),
    variant_id_b38 != ""
  ) %>%
  distinct()

cbpnet_vars <- cbpnet_scores_clean %>%
  inner_join(
    cbpnet_mapping_clean,
    by = "variant_id"
  )

if (nrow(cbpnet_vars) == 0L) {
  stop(
    "No ChromBPNet score variants matched the mapping table by variant_id."
  )
}

unmapped_cbpnet_scores <- cbpnet_scores_clean %>%
  anti_join(
    cbpnet_mapping_clean %>% distinct(variant_id),
    by = "variant_id"
  )

write_csv(
  unmapped_cbpnet_scores,
  file.path(
    output_dir,
    "QC_ChromBPNet_score_variants_without_hg19_mapping.csv"
  )
)

message(
  "ChromBPNet score variants: ",
  nrow(cbpnet_scores_clean)
)

message(
  "ChromBPNet score variants mapped to hg19: ",
  n_distinct(cbpnet_vars$variant_id)
)

# ======================================================================
# 3) LOAD QTL CREDIBLE SETS
# ======================================================================

eqtl_cs <- load_cs_files(
  eqtl_cs_dir,
  "^All_CD4T_cells_chr[0-9]+_credible_sets\\.txt$"
)

caqtl_cs <- load_cs_files(
  caqtl_cs_dir,
  "^CD4T_chromatin_chr[0-9]+_credible_sets\\.txt$"
)

required_eqtl_cols <- c(
  "gene",
  "region",
  "cs_index",
  "variant_id"
)

required_caqtl_cols <- c(
  "peak",
  "region",
  "cs_index",
  "variant_id"
)

missing_eqtl_cols <- setdiff(
  required_eqtl_cols,
  names(eqtl_cs)
)

missing_caqtl_cols <- setdiff(
  required_caqtl_cols,
  names(caqtl_cs)
)

if (length(missing_eqtl_cols) > 0L) {
  stop(
    "eQTL CS data missing: ",
    paste(missing_eqtl_cols, collapse = ", ")
  )
}

if (length(missing_caqtl_cols) > 0L) {
  stop(
    "caQTL CS data missing: ",
    paste(missing_caqtl_cols, collapse = ", ")
  )
}

qtl_cs_lookup <- bind_rows(
  eqtl_cs %>%
    group_by(gene, region, cs_index) %>%
    summarise(
      n_qtl_variants = n(),
      .groups = "drop"
    ) %>%
    transmute(
      qtl_type = "eQTL",
      feature = as.character(gene),
      qtl_region = as.character(region),
      qtl_cs_index = as.integer(cs_index),
      n_qtl_variants
    ),
  
  caqtl_cs %>%
    group_by(peak, region, cs_index) %>%
    summarise(
      n_qtl_variants = n(),
      .groups = "drop"
    ) %>%
    transmute(
      qtl_type = "caQTL",
      feature = as.character(peak),
      qtl_region = as.character(region),
      qtl_cs_index = as.integer(cs_index),
      n_qtl_variants
    )
)

# ======================================================================
# 4) LOAD GWAS FINE-MAPPED CREDIBLE SETS
# ======================================================================
#
# This preserves actual SuSiE component numbers from names such as
# L1, L2, L5, and L9. It does not replace them with list positions.

extract_gwas_cs_variants <- function(rds_path) {
  region_obj <- readRDS(rds_path)
  
  if (isFALSE(region_obj$converged)) {
    return(NULL)
  }
  
  cs_list <- region_obj$sets$cs
  
  if (is.null(cs_list) || length(cs_list) == 0L) {
    return(NULL)
  }
  
  all_vars <- colnames(region_obj$alpha)
  
  if (is.null(all_vars)) {
    stop(
      "No variant column names found in alpha matrix: ",
      rds_path
    )
  }
  
  cs_names <- names(cs_list)
  
  if (
    is.null(cs_names) ||
    any(is.na(cs_names)) ||
    any(cs_names == "")
  ) {
    stop(
      "SuSiE credible sets lack component names in: ",
      rds_path
    )
  }
  
  cs_indices <- suppressWarnings(
    as.integer(sub("^L", "", cs_names))
  )
  
  if (anyNA(cs_indices)) {
    stop(
      "Could not parse SuSiE component indices from: ",
      paste(cs_names, collapse = ", "),
      "\nFile: ",
      rds_path
    )
  }
  
  bind_rows(
    lapply(seq_along(cs_list), function(i) {
      idx <- cs_list[[i]]
      
      tibble(
        variant_id = all_vars[idx],
        cs_index = cs_indices[i],
        pip = region_obj$pip[idx]
      )
    })
  )
}

load_gwas_cs_variants <- function(gwas_id) {
  message("Loading GWAS CSs: ", gwas_id)
  
  files <- list.files(
    path = file.path(
      gwas_fm_dir,
      "output/nathan_completed",
      gwas_id
    ),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  files <- files[
    grepl("_1\\.5Mb/susie/", files)
  ]
  
  if (length(files) == 0L) {
    return(NULL)
  }
  
  out <- bind_rows(
    lapply(files, function(f) {
      df <- extract_gwas_cs_variants(f)
      
      if (is.null(df)) {
        return(NULL)
      }
      
      df$region <- paste0(
        "chr",
        sub(
          ".*_([0-9]+\\.[0-9]+)\\.rds$",
          "\\1",
          f
        )
      )
      
      df
    })
  )
  
  if (nrow(out) == 0L) {
    return(NULL)
  }
  
  out %>%
    mutate(
      gwas_id = gwas_id,
      trait = extract_trait(gwas_id),
      region_cs = make_region_cs(
        region,
        cs_index
      )
    )
}

file_names <- unlist(
  fread(
    file.path(
      data_dir,
      "preprocessed_folder_names.txt"
    ),
    header = FALSE
  ),
  use.names = FALSE
)

file_names <- unique(
  as.character(file_names)
)

gwas_variant_all <- bind_rows(
  lapply(
    file_names,
    load_gwas_cs_variants
  )
)

if (nrow(gwas_variant_all) == 0L) {
  stop("No GWAS credible sets were loaded.")
}

gwas_cs_study <- gwas_variant_all %>%
  group_by(
    gwas_id,
    trait,
    region,
    cs_index,
    region_cs
  ) %>%
  summarise(
    n_gwas_variants = n_distinct(variant_id),
    .groups = "drop"
  )

# ======================================================================
# 5) LOAD SIGNIFICANT GWAS-QTL COLOCALIZATION PAIRS
# ======================================================================

coloc_eqtl_dir <- file.path(
  data_dir,
  "coloc_results/eqtl_gwas_coloc"
)

coloc_caqtl_dir <- file.path(
  data_dir,
  "coloc_results"
)

read_coloc_pairs <- function(
    path,
    gwas_id,
    qtl_type,
    feature_col,
    qtl_cs,
    qtl_region_col = "region_GRC38"
) {
  if (!file.exists(path)) {
    return(NULL)
  }
  
  dt <- fread(
    path,
    showProgress = FALSE
  )
  
  if (nrow(dt) == 0L) {
    return(NULL)
  }
  
  needed <- c(
    "pval",
    "pval_nominal",
    "region",
    "idx1",
    "idx2",
    feature_col,
    qtl_region_col
  )
  
  missing <- setdiff(
    needed,
    names(dt)
  )
  
  if (length(missing) > 0L) {
    stop(
      qtl_type,
      " coloc file missing: ",
      paste(missing, collapse = ", "),
      "\nFile: ",
      path
    )
  }
  
  pp4_candidates <- c(
    "PP.H4.abf",
    "PP.H4",
    "PP.H4.abf.x",
    "PP.H4.abf.y"
  )
  
  pp4_col <- pp4_candidates[
    pp4_candidates %in% names(dt)
  ][1]
  
  if (length(pp4_col) == 0L || is.na(pp4_col)) {
    stop(
      "No PP.H4 column found in coloc file: ",
      path,
      ". A PP.H4 threshold is required for this analysis."
    )
  }
  
  dt[, coloc_PP.H4.abf := as.numeric(get(pp4_col))]
  
  dt <- dt[
    pval < GWAS_P_THRESH &
      pval_nominal < QTL_P_THRESH &
      coloc_PP.H4.abf > PP_H4_THRESH
  ]
  
  if (nrow(dt) == 0L) {
    return(NULL)
  }
  
  out <- dt %>%
    as_tibble() %>%
    transmute(
      gwas_id = gwas_id,
      trait = extract_trait(gwas_id),
      
      region = as.character(region),
      cs_index = as.integer(idx1),
      region_cs = make_region_cs(
        region,
        idx1
      ),
      
      qtl_type = qtl_type,
      feature = as.character(
        .data[[feature_col]]
      ),
      qtl_region = as.character(
        .data[[qtl_region_col]]
      ),
      qtl_cs_index = as.integer(idx2),
      
      coloc_PP.H4.abf = coloc_PP.H4.abf,
      coloc_source_file = path
    ) %>%
    distinct()
  
  qtl_feature_lookup_col <- if (
    qtl_type == "eQTL"
  ) {
    "gene"
  } else {
    "peak"
  }
  
  qtl_key_check <- out %>%
    distinct(
      qtl_type,
      feature,
      qtl_region,
      qtl_cs_index
    ) %>%
    left_join(
      qtl_cs %>%
        transmute(
          feature = as.character(
            .data[[qtl_feature_lookup_col]]
          ),
          qtl_region = as.character(region),
          qtl_cs_index = as.integer(cs_index),
          qtl_cs_lookup_matched = TRUE
        ) %>%
        distinct(),
      by = c(
        "feature",
        "qtl_region",
        "qtl_cs_index"
      )
    ) %>%
    mutate(
      qtl_cs_lookup_matched = replace_na(
        qtl_cs_lookup_matched,
        FALSE
      )
    )
  
  if (any(!qtl_key_check$qtl_cs_lookup_matched)) {
    warning(
      qtl_type,
      " coloc rows from ",
      basename(path),
      " include QTL CS keys not found in the loaded CS table."
    )
  }
  
  out
}

load_study_coloc_pairs <- function(gwas_id) {
  eqtl_path <- file.path(
    coloc_eqtl_dir,
    paste0(
      gwas_id,
      "_coloc_results.csv"
    )
  )
  
  caqtl_path <- file.path(
    coloc_caqtl_dir,
    paste0(
      "CD4T_chromatin_",
      gwas_id,
      "_coloc_results.csv"
    )
  )
  
  bind_rows(
    read_coloc_pairs(
      eqtl_path,
      gwas_id,
      "eQTL",
      "gene",
      eqtl_cs,
      qtl_region_col = "region_GRC38"
    ),
    
    read_coloc_pairs(
      caqtl_path,
      gwas_id,
      "caQTL",
      "peak",
      caqtl_cs,
      qtl_region_col = "region_GRC38"
    )
  )
}

coloc_pairs_raw <- bind_rows(
  lapply(
    file_names,
    load_study_coloc_pairs
  )
) %>%
  distinct()

if (nrow(coloc_pairs_raw) == 0L) {
  warning(
    "No significant coloc pairs were loaded. ",
    "All loci will be classified as no_coloc."
  )
}

# ======================================================================
# 6) COLOC AND QTL-CS QC
# ======================================================================

coloc_qtl_key_qc <- coloc_pairs_raw %>%
  distinct(
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index
  ) %>%
  left_join(
    qtl_cs_lookup %>%
      select(
        all_of(KEY_QTL),
        n_qtl_variants
      ) %>%
      distinct(),
    by = KEY_QTL
  ) %>%
  mutate(
    qtl_cs_key_matched = !is.na(
      n_qtl_variants
    )
  )

print(
  coloc_qtl_key_qc %>%
    count(
      qtl_type,
      qtl_cs_key_matched
    )
)

unmatched_coloc_qtl_keys <- coloc_qtl_key_qc %>%
  filter(!qtl_cs_key_matched)

write_csv(
  unmatched_coloc_qtl_keys,
  file.path(
    output_dir,
    "QC_unmatched_coloc_qtl_keys_after_parsing.csv"
  )
)

# ======================================================================
# 7) REQUIRE STUDY-SPECIFIC GWAS-CS MATCHING
# ======================================================================

gwas_study_key <- gwas_cs_study %>%
  select(all_of(KEY_STUDY)) %>%
  distinct()

excluded_coloc_pairs_not_in_primary_gwas_cs <- coloc_pairs_raw %>%
  anti_join(
    gwas_study_key,
    by = KEY_STUDY
  ) %>%
  arrange(
    trait,
    region,
    cs_index,
    gwas_id,
    qtl_type,
    feature
  )

write_csv(
  excluded_coloc_pairs_not_in_primary_gwas_cs,
  file.path(
    output_dir,
    "QC_coloc_pairs_without_matching_study_specific_GWAS_CS.csv"
  )
)

message(
  "Significant coloc rows excluded because their study-specific GWAS CS ",
  "was absent: ",
  nrow(excluded_coloc_pairs_not_in_primary_gwas_cs)
)

coloc_pairs_primary_study <- coloc_pairs_raw %>%
  inner_join(
    gwas_study_key,
    by = KEY_STUDY
  ) %>%
  left_join(
    qtl_cs_lookup,
    by = KEY_QTL
  ) %>%
  mutate(
    qtl_cs_matched = !is.na(
      n_qtl_variants
    )
  )

unmatched_primary_pairs_to_qtl_cs <- coloc_pairs_primary_study %>%
  filter(!qtl_cs_matched)

write_csv(
  unmatched_primary_pairs_to_qtl_cs,
  file.path(
    output_dir,
    "QC_primary_coloc_pairs_without_QTL_CS.csv"
  )
)

message(
  "Primary coloc rows without a matched QTL CS: ",
  nrow(unmatched_primary_pairs_to_qtl_cs)
)

# Use only coloc pairs whose GWAS CS and QTL CS keys were validated.
coloc_pairs_primary_study_valid <- coloc_pairs_primary_study %>%
  filter(qtl_cs_matched)

# ======================================================================
# 8) DEFINE ONE CANONICAL TRAIT-GWAS LOCUS RECORD
# ======================================================================
#
# NOTE:
# This retains your existing trait + region + cs_index canonicalization.
# If the same trait-region-CS index appears across multiple studies, the
# studies are collapsed. The QC table below records these repeated loci.

coloc_cs_study <- coloc_pairs_primary_study_valid %>%
  group_by(
    gwas_id,
    trait,
    region,
    cs_index,
    region_cs
  ) %>%
  summarise(
    eqtl_coloc = any(qtl_type == "eQTL"),
    caqtl_coloc = any(qtl_type == "caQTL"),
    .groups = "drop"
  )

gwas_cs_study_annot <- gwas_cs_study %>%
  left_join(
    coloc_cs_study,
    by = KEY_STUDY
  ) %>%
  mutate(
    eqtl_coloc = replace_na(
      eqtl_coloc,
      FALSE
    ),
    caqtl_coloc = replace_na(
      caqtl_coloc,
      FALSE
    )
  )

trait_locus_master <- gwas_cs_study_annot %>%
  group_by(
    trait,
    region,
    cs_index,
    region_cs
  ) %>%
  summarise(
    n_studies_at_locus = n_distinct(gwas_id),
    study_ids = paste(
      sort(unique(gwas_id)),
      collapse = ";"
    ),
    eqtl_coloc = any(eqtl_coloc),
    caqtl_coloc = any(caqtl_coloc),
    .groups = "drop"
  ) %>%
  mutate(
    qtl_coloc_category = case_when(
      eqtl_coloc & caqtl_coloc ~ "both",
      !eqtl_coloc & caqtl_coloc ~ "just_caqtl",
      eqtl_coloc & !caqtl_coloc ~ "just_eqtl",
      TRUE ~ "no_coloc"
    )
  )

repeated_trait_loci <- trait_locus_master %>%
  filter(n_studies_at_locus > 1L)

write_csv(
  repeated_trait_loci,
  file.path(
    output_dir,
    "QC_repeated_trait_loci_collapsed_across_studies.csv"
  )
)

message(
  "Repeated trait-loci collapsed across studies: ",
  nrow(repeated_trait_loci)
)

stopifnot(
  nrow(
    trait_locus_master %>%
      count(
        trait,
        region,
        cs_index,
        region_cs
      ) %>%
      filter(n > 1L)
  ) == 0L
)

# ======================================================================
# 9) BUILD UNAMBIGUOUS CHROMBPNET LOOKUP
# ======================================================================
#
# Exclude hg19 keys mapping to multiple hg38 variants in the primary test.

cbpnet_lookup_all <- cbpnet_vars %>%
  transmute(
    gwas_variant_hg19 = normalize_hg19_key(
      hg19_var_lookup
    ),
    variant_key_hg38 = normalize_hg38_key(
      variant_id
    ),
    chrombpnet_hit = replace_na(
      as.logical(chrombpnet_hit),
      FALSE
    )
  ) %>%
  filter(
    !is.na(gwas_variant_hg19),
    gwas_variant_hg19 != "",
    !is.na(variant_key_hg38),
    variant_key_hg38 != ""
  ) %>%
  group_by(
    gwas_variant_hg19,
    variant_key_hg38
  ) %>%
  summarise(
    chrombpnet_hit = any(
      chrombpnet_hit
    ),
    .groups = "drop"
  )

ambiguous_hg19_keys <- cbpnet_lookup_all %>%
  distinct(
    gwas_variant_hg19,
    variant_key_hg38
  ) %>%
  count(
    gwas_variant_hg19,
    name = "n_hg38_mappings"
  ) %>%
  filter(
    n_hg38_mappings > 1L
  )

write_csv(
  ambiguous_hg19_keys,
  file.path(
    output_dir,
    "QC_ambiguous_hg19_to_hg38_keys_excluded.csv"
  )
)

cbpnet_lookup_primary <- cbpnet_lookup_all %>%
  anti_join(
    ambiguous_hg19_keys %>%
      select(gwas_variant_hg19),
    by = "gwas_variant_hg19"
  ) %>%
  mutate(
    chrombpnet_callable = TRUE
  )

message(
  "Ambiguous hg19 keys excluded: ",
  nrow(ambiguous_hg19_keys)
)

# ======================================================================
# 10) COUNT CALLABLE/HIT VARIANTS IN GWAS CREDIBLE SETS ONLY
# ======================================================================

gwas_cs_callable_variants <- gwas_variant_all %>%
  transmute(
    trait = as.character(trait),
    region = as.character(region),
    cs_index = as.integer(cs_index),
    region_cs = as.character(region_cs),
    gwas_id = as.character(gwas_id),
    gwas_variant_hg19 = normalize_hg19_key(
      variant_id
    ),
    variant_chr = parse_variant_chr(
      variant_id
    ),
    variant_pos = parse_variant_pos(
      variant_id
    )
  ) %>%
  left_join(
    cbpnet_lookup_primary,
    by = "gwas_variant_hg19",
    relationship = "many-to-one"
  ) %>%
  filter(
    !is.na(variant_key_hg38)
  ) %>%
  distinct(
    across(all_of(KEY_LOCUS)),
    variant_key_hg38,
    .keep_all = TRUE
  )

# Derive MHC overlap from actual GWAS-CS variant positions.
gwas_locus_bounds <- gwas_variant_all %>%
  transmute(
    trait = as.character(trait),
    region = as.character(region),
    cs_index = as.integer(cs_index),
    region_cs = as.character(region_cs),
    variant_chr = parse_variant_chr(
      variant_id
    ),
    variant_pos = parse_variant_pos(
      variant_id
    )
  ) %>%
  filter(
    !is.na(variant_chr),
    !is.na(variant_pos)
  ) %>%
  group_by(
    across(all_of(KEY_LOCUS))
  ) %>%
  summarise(
    chromosome = first(variant_chr),
    locus_start = min(
      variant_pos,
      na.rm = TRUE
    ),
    locus_end = max(
      variant_pos,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    extended_mhc = (
      chromosome == MHC_CHR &
        locus_start <= MHC_END &
        locus_end >= MHC_START
    )
  )

gwas_only_counts <- gwas_cs_callable_variants %>%
  group_by(
    across(all_of(KEY_LOCUS))
  ) %>%
  summarise(
    m_callable = n_distinct(
      variant_key_hg38
    ),
    
    h_hits = n_distinct(
      variant_key_hg38[
        chrombpnet_hit %in% TRUE
      ]
    ),
    
    callable_variants = collapse_unique_values(
      variant_key_hg38
    ),
    
    prioritized_variants = collapse_unique_values(
      variant_key_hg38[
        chrombpnet_hit %in% TRUE
      ]
    ),
    
    .groups = "drop"
  )

analysis_master <- trait_locus_master %>%
  left_join(
    gwas_only_counts,
    by = KEY_LOCUS
  ) %>%
  left_join(
    gwas_locus_bounds,
    by = KEY_LOCUS
  ) %>%
  mutate(
    m_callable = replace_na(
      as.integer(m_callable),
      0L
    ),
    
    h_hits = replace_na(
      as.integer(h_hits),
      0L
    ),
    
    frac_hits = if_else(
      m_callable > 0L,
      h_hits / m_callable,
      NA_real_
    ),
    
    callable_bin = make_callable_bin(
      m_callable
    ),
    
    extended_mhc = replace_na(
      extended_mhc,
      FALSE
    )
  )

stopifnot(
  all(
    analysis_master$h_hits <=
      analysis_master$m_callable
  ),
  
  nrow(
    analysis_master %>%
      count(
        across(all_of(KEY_LOCUS))
      ) %>%
      filter(n > 1L)
  ) == 0L
)

write_csv(
  analysis_master,
  file.path(
    output_dir,
    "GWAS_CS_only_callable_hit_counts_by_locus.csv"
  )
)

analysis_loci <- analysis_master %>%
  filter(
    m_callable >= 1L,
    qtl_coloc_category %in% c(
      "no_coloc",
      "both",
      "just_caqtl",
      "just_eqtl"
    ),
    !is.na(chromosome),
    !is.na(callable_bin)
  )

if (EXCLUDE_EXTENDED_MHC) {
  message(
    "Callable loci excluded for extended MHC overlap: ",
    sum(analysis_loci$extended_mhc)
  )
  
  analysis_loci <- analysis_loci %>%
    filter(!extended_mhc)
}

analysis_loci <- analysis_loci %>%
  mutate(
    qtl_coloc_category = factor(
      qtl_coloc_category,
      levels = c(
        "no_coloc",
        "both",
        "just_caqtl",
        "just_eqtl"
      )
    ),
    
    chromosome = factor(
      chromosome
    ),
    
    callable_bin = droplevels(
      callable_bin
    )
  )

message(
  "Callable canonical loci retained: ",
  nrow(analysis_loci)
)

category_summary <- analysis_loci %>%
  group_by(
    qtl_coloc_category
  ) %>%
  summarise(
    n_loci = n(),
    total_callable = sum(
      m_callable
    ),
    total_hits = sum(
      h_hits
    ),
    pooled_hit_fraction = safe_fraction(
      h_hits,
      m_callable
    ),
    mean_locus_hit_fraction = mean(
      frac_hits,
      na.rm = TRUE
    ),
    median_callable_per_locus = median(
      m_callable
    ),
    .groups = "drop"
  )

print(category_summary)

write_csv(
  category_summary,
  file.path(
    output_dir,
    "GWAS_CS_only_category_summary.csv"
  )
)

# ======================================================================
# 11) PRIMARY LOCUS-AWARE QUASI-BINOMIAL MODEL
# ======================================================================
#
# One row per canonical GWAS locus.
#
# Response:
#   h_hits prioritized variants
#   m_callable - h_hits non-prioritized callable variants
#
# Predictors:
#   qtl_coloc_category
#   chromosome
#   callable count stratum

quasi_fit <- glm(
  cbind(
    h_hits,
    m_callable - h_hits
  ) ~
    qtl_coloc_category +
    chromosome +
    callable_bin,
  family = quasibinomial(),
  data = analysis_loci
)

saveRDS(
  quasi_fit,
  file.path(
    output_dir,
    "GWAS_CS_only_quasibinomial_model.rds"
  )
)

coef_table <- summary(
  quasi_fit
)$coefficients

extract_quasi_contrast <- function(category_label) {
  coef_name <- paste0(
    "qtl_coloc_category",
    category_label
  )
  
  if (!coef_name %in% rownames(coef_table)) {
    stop(
      "Coefficient not estimable: ",
      coef_name
    )
  }
  
  estimate <- unname(
    coef_table[
      coef_name,
      "Estimate"
    ]
  )
  
  se <- unname(
    coef_table[
      coef_name,
      "Std. Error"
    ]
  )
  
  p <- unname(
    coef_table[
      coef_name,
      "Pr(>|t|)"
    ]
  )
  
  tibble(
    category = category_label,
    model_logOR = estimate,
    model_OR = exp(estimate),
    model_CI_low = exp(
      estimate - 1.96 * se
    ),
    model_CI_high = exp(
      estimate + 1.96 * se
    ),
    model_p_two_sided = p
  )
}

model_results <- bind_rows(
  extract_quasi_contrast(
    "both"
  ),
  extract_quasi_contrast(
    "just_caqtl"
  ),
  extract_quasi_contrast(
    "just_eqtl"
  )
)

# ======================================================================
# 12) STRATIFIED PERMUTATION TEST USING POOLED LOG ODDS RATIO
# ======================================================================

safe_log_or <- function(
    category,
    h_hits,
    m_callable,
    target_category,
    baseline_category = "no_coloc",
    correction = 0.5
) {
  category <- as.character(
    category
  )
  
  target <- (
    category == target_category
  )
  
  baseline <- (
    category == baseline_category
  )
  
  if (!any(target) || !any(baseline)) {
    return(NA_real_)
  }
  
  target_hits <- sum(
    h_hits[target],
    na.rm = TRUE
  )
  
  target_nonhits <- sum(
    m_callable[target] -
      h_hits[target],
    na.rm = TRUE
  )
  
  baseline_hits <- sum(
    h_hits[baseline],
    na.rm = TRUE
  )
  
  baseline_nonhits <- sum(
    m_callable[baseline] -
      h_hits[baseline],
    na.rm = TRUE
  )
  
  log(
    (
      (target_hits + correction) /
        (target_nonhits + correction)
    ) /
      (
        (baseline_hits + correction) /
          (baseline_nonhits + correction)
      )
  )
}

permute_labels_within_strata <- function(
    category,
    chromosome,
    callable_bin
) {
  category <- as.character(
    category
  )
  
  strata <- interaction(
    chromosome,
    callable_bin,
    drop = TRUE,
    lex.order = TRUE
  )
  
  out <- category
  
  for (s in levels(strata)) {
    idx <- which(
      strata == s
    )
    
    if (
      length(idx) > 1L &&
      length(
        unique(category[idx])
      ) > 1L
    ) {
      out[idx] <- sample(
        category[idx],
        replace = FALSE
      )
    }
  }
  
  out
}

run_pairwise_permutation <- function(
    locus_df,
    target_category,
    nperm = 100000L,
    seed = 1L
) {
  dat <- locus_df %>%
    filter(
      qtl_coloc_category %in%
        c(
          "no_coloc",
          target_category
        )
    ) %>%
    mutate(
      category_chr = as.character(
        qtl_coloc_category
      ),
      chromosome_chr = as.character(
        chromosome
      )
    )
  
  target <- (
    dat$category_chr ==
      target_category
  )
  
  baseline <- (
    dat$category_chr ==
      "no_coloc"
  )
  
  if (!any(target) || !any(baseline)) {
    stop(
      "Missing target or no_coloc loci for: ",
      target_category
    )
  }
  
  observed_logOR <- safe_log_or(
    category = dat$category_chr,
    h_hits = dat$h_hits,
    m_callable = dat$m_callable,
    target_category = target_category
  )
  
  observed_OR <- exp(
    observed_logOR
  )
  
  stratum_qc <- dat %>%
    mutate(
      stratum = interaction(
        chromosome_chr,
        callable_bin,
        drop = TRUE
      )
    ) %>%
    group_by(stratum) %>%
    summarise(
      n_loci = n(),
      n_categories = n_distinct(
        category_chr
      ),
      .groups = "drop"
    )
  
  n_movable_loci <- sum(
    stratum_qc$n_loci[
      stratum_qc$n_categories > 1L
    ]
  )
  
  if (n_movable_loci == 0L) {
    stop(
      "No labels can move within chromosome + callable-bin strata for ",
      target_category,
      ". Broaden CALLABLE_BREAKS."
    )
  }
  
  set.seed(seed)
  
  null_logOR <- replicate(
    nperm,
    {
      permuted_category <- permute_labels_within_strata(
        category = dat$category_chr,
        chromosome = dat$chromosome_chr,
        callable_bin = dat$callable_bin
      )
      
      safe_log_or(
        category = permuted_category,
        h_hits = dat$h_hits,
        m_callable = dat$m_callable,
        target_category = target_category
      )
    }
  )
  
  null_logOR <- null_logOR[
    is.finite(null_logOR)
  ]
  
  if (length(null_logOR) == 0L) {
    stop(
      "No finite permutation statistics for: ",
      target_category
    )
  }
  
  p_enrichment <- (
    sum(
      null_logOR >= observed_logOR
    ) + 1
  ) / (
    length(null_logOR) + 1
  )
  
  p_depletion <- (
    sum(
      null_logOR <= observed_logOR
    ) + 1
  ) / (
    length(null_logOR) + 1
  )
  
  p_two_sided <- (
    sum(
      abs(null_logOR) >=
        abs(observed_logOR)
    ) + 1
  ) / (
    length(null_logOR) + 1
  )
  
  tibble(
    category = target_category,
    
    n_target_loci = sum(target),
    n_baseline_loci = sum(baseline),
    n_movable_loci = n_movable_loci,
    
    target_total_callable = sum(
      dat$m_callable[target]
    ),
    
    target_total_hits = sum(
      dat$h_hits[target]
    ),
    
    baseline_total_callable = sum(
      dat$m_callable[baseline]
    ),
    
    baseline_total_hits = sum(
      dat$h_hits[baseline]
    ),
    
    target_pooled_hit_fraction = safe_fraction(
      dat$h_hits[target],
      dat$m_callable[target]
    ),
    
    baseline_pooled_hit_fraction = safe_fraction(
      dat$h_hits[baseline],
      dat$m_callable[baseline]
    ),
    
    observed_pooled_logOR = observed_logOR,
    observed_pooled_OR = observed_OR,
    
    permutation_null_mean_logOR = mean(
      null_logOR
    ),
    
    permutation_null_sd_logOR = sd(
      null_logOR
    ),
    
    perm_p_enrichment = p_enrichment,
    perm_p_depletion = p_depletion,
    perm_p_two_sided = p_two_sided,
    
    nperm = length(null_logOR),
    seed = seed,
    
    null_logOR = list(
      null_logOR
    )
  )
}

permutation_results <- bind_rows(
  run_pairwise_permutation(
    analysis_loci,
    "both",
    N_PERM,
    PERM_SEEDS[["both"]]
  ),
  
  run_pairwise_permutation(
    analysis_loci,
    "just_caqtl",
    N_PERM,
    PERM_SEEDS[["just_caqtl"]]
  ),
  
  run_pairwise_permutation(
    analysis_loci,
    "just_eqtl",
    N_PERM,
    PERM_SEEDS[["just_eqtl"]]
  )
)

saveRDS(
  permutation_results,
  file.path(
    output_dir,
    "GWAS_CS_only_permutation_results_with_nulls.rds"
  )
)

# ======================================================================
# 13) COMBINE MODEL AND PERMUTATION RESULTS
# ======================================================================

final_results <- model_results %>%
  left_join(
    permutation_results %>%
      select(-null_logOR),
    by = "category"
  ) %>%
  mutate(
    observed_direction = case_when(
      model_OR > 1 ~ "enrichment",
      model_OR < 1 ~ "depletion",
      TRUE ~ "null"
    ),
    
    perm_p_directional = case_when(
      observed_direction == "enrichment" ~
        perm_p_enrichment,
      
      observed_direction == "depletion" ~
        perm_p_depletion,
      
      TRUE ~
        perm_p_two_sided
    )
  )

print(final_results)

write_csv(
  final_results,
  file.path(
    output_dir,
    "GWAS_CS_only_correct_enrichment_results.csv"
  )
)

# ======================================================================
# 14) GLOBAL TEST OF THE FOUR-LEVEL COLOCALIZATION FACTOR
# ======================================================================

quasi_fit_null <- update(
  quasi_fit,
  . ~ . - qtl_coloc_category
)

global_test <- anova(
  quasi_fit_null,
  quasi_fit,
  test = "F"
)

capture.output(
  global_test,
  file = file.path(
    output_dir,
    "GWAS_CS_only_global_colocalization_factor_test.txt"
  )
)

print(global_test)

# ======================================================================
# 15) FOREST PLOT
# ======================================================================

plot_df <- final_results %>%
  mutate(
    category_label = recode(
      category,
      both = "eQTL + caQTL",
      just_caqtl = "caQTL only",
      just_eqtl = "eQTL only"
    ),
    
    category_label = factor(
      category_label,
      levels = rev(
        c(
          "eQTL + caQTL",
          "caQTL only",
          "eQTL only"
        )
      )
    )
  )

forest_plot <- ggplot(
  plot_df,
  aes(
    x = model_OR,
    y = category_label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = model_CI_low,
      xmax = model_CI_high
    ),
    height = 0,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_log10(
    breaks = log_breaks(
      n = 6
    ),
    labels = label_number(
      accuracy = 0.01
    )
  ) +
  labs(
    x = paste0(
      "Odds ratio for ChromBPNet prioritization\n",
      "among callable variants in GWAS credible sets"
    ),
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  )

print(forest_plot)

ggsave(
  filename = file.path(
    output_dir,
    "GWAS_CS_only_correct_enrichment_forest_plot.pdf"
  ),
  plot = forest_plot,
  width = 7.5,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = file.path(
    output_dir,
    "GWAS_CS_only_correct_enrichment_forest_plot.png"
  ),
  plot = forest_plot,
  width = 7.5,
  height = 4.5,
  units = "in",
  dpi = 400
)

message("Analysis completed.")
message(
  "Primary results: ",
  file.path(
    output_dir,
    "GWAS_CS_only_correct_enrichment_results.csv"
  )
)


