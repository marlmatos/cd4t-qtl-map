## ======================================================================

## ChromBPNet support at GWAS loci and QTL-GWAS colocalized loci

## Unified workflow for plots + Supplementary Table 8

## ======================================================================

##

## Analytical unit:

## One unique GWAS fine-mapped credible-set locus per trait:

## trait + region + cs_index + region_cs

##

## Canonical algorithm:

## 1. Load GWAS credible sets and define the study-specific GWAS-CS universe.

## 2. Retain significant coloc pairs only when the matching GWAS CS exists

## in the same GWAS study.

## 3. Annotate ChromBPNet support using the same validated membership calls

## throughout the workflow:

## - GWAS CS: variant_id %in% cbpnet_set37

## - QTL  CS: variant_id %in% cbpnet_set_b38

## 4. Build `trait_locus_master`, the canonical object used to classify loci

## for the plots and determine Supplementary Table 8 inclusion.

## 5. Build plots by summarising `trait_locus_master`, ensuring that each

## trait-GWAS locus contributes only once.

## 6. Select the exact QTL-colocalized loci with ChromBPNet support from

## `trait_locus_master`.

## 7. Collapse all linked eGenes, caPeaks, and paired QTL credible sets at

## each retained locus into a single locus-level row.

## 8. Export Supplementary Table 8 with exactly one row per unique

## trait-GWAS locus, matching the analytical unit used in the plots.

##

## Interpretation of Supplementary Table 8:

## A retained locus has QTL-GWAS colocalization and ChromBPNet support in

## its GWAS credible set, at least one paired QTL credible set, or both.

##

## Each retained trait-GWAS locus is reported once. When multiple eQTL or

## caQTL credible sets colocalize with the same GWAS locus, their features

## and credible-set identifiers are collapsed into semicolon-separated

## fields rather than reported as duplicated rows.

##

## The table distinguishes:

## - all paired QTL credible sets at the locus;

## - the subset of paired QTL credible sets that themselves contain

## a ChromBPNet-supported variant;

## - ChromBPNet-supported variants in the GWAS credible set;

## - ChromBPNet-supported variants in supporting paired QTL credible sets.

##

## Locus-level counts in the plots and Supplementary Table 8 must therefore

## agree exactly, and duplicate trait-GWAS locus rows are treated as an error.

##

## Edit only the CONFIGURATION section before running.

## ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

## ======================================================================
## 0) CONFIGURATION: EDIT THESE PATHS
## ======================================================================
source("~/cd4t-qtl-map/analysis/utils/color_pallete_helper.R")
source("~/cd4t-qtl-map/analysis/paths_config.R")

## Base data directory containing coloc_results/ and preprocessed_folder_names.txt
data.dir       <- script12$data.dir

## ChromBPNet-prioritized variant table.
## Required columns used below:
##   hg19_chr, hg19_pos, alt, ref, variant_id
## Optional display column: rsid
chrombpnet_file <- script12$chrombpnet_file

## Molecular QTL credible set directories.
eqtl_cs_dir    <- script12$eqtl_cs_dir
caqtl_cs_dir   <- script12$caqtl_cs_dir

## GWAS fine-mapping root directory containing output/nathan_completed/<gwas_id>/...
gwas_fm_dir    <- script12$gwas_fm_dir  # SuSiE .rds files

## Outputs
output_dir <-  script12$output_dir
plot_output_file <- script12$plot_output_file
output_file    <- script12$output_file


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_output_file), recursive = TRUE, showWarnings = FALSE)

## Coloc significance thresholds retained from the original plot workflow.
GWAS_P_THRESH <- 1e-5
QTL_P_THRESH  <- 1e-3
PP_H4_THRESH   <- 0.5
PIP_THRESH     <- 0.0  # filter QTL/GWAS CS variants by PIP
## Optional pretty trait labels. Add mappings here if desired.
## Example: trait_labels <- c("IBD" = "Inflammatory bowel disease")
trait_labels <- character(0)

## Fail immediately if configuration paths have not been edited.
configured_paths <- c(data.dir, chrombpnet_file, eqtl_cs_dir, caqtl_cs_dir)
if (any(grepl("REPLACE_WITH", configured_paths, fixed = TRUE))) {
  stop(
    "Edit the CONFIGURATION section: set data.dir, chrombpnet_file, ",
    "eqtl_cs_dir, and caqtl_cs_dir before running this script."
  )
}

required_files <- c(
  chrombpnet_file,
  file.path(data.dir, "preprocessed_folder_names.txt")
)
required_dirs <- c(
  data.dir,
  eqtl_cs_dir,
  caqtl_cs_dir,
  gwas_fm_dir
)
if (any(!file.exists(required_files))) {
  stop("Missing input file(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}
if (any(!dir.exists(required_dirs))) {
  stop("Missing input directory/directories: ", paste(required_dirs[!dir.exists(required_dirs)], collapse = ", "))
}
## ======================================================================
## 1) CONSTANT JOIN KEYS AND HELPERS
## ======================================================================

KEY_STUDY <- c("gwas_id", "trait", "region", "cs_index", "region_cs")
KEY_LOCUS <- c("trait", "region", "cs_index", "region_cs")
KEY_QTL   <- c("qtl_type", "feature", "qtl_region", "qtl_cs_index")

source("~/cd4t-qtl-map/analysis/utils/color_pallete_helper.R")

label_trait <- function(x) dplyr::recode(x, !!!trait_labels, .default = x)

extract_trait <- function(gwas_id) {
  x <- sub("_preprocessed$", "", gwas_id)
  pieces <- strsplit(x, "_", fixed = TRUE)[[1]]
  if (length(pieces) < 3) stop("Cannot parse trait from GWAS ID: ", gwas_id)
  pieces[3]
}

make_region_cs <- function(region, cs_index) paste(region, cs_index, sep = "_")

collapse_unique_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  x <- sort(unique(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

collapse_variant_strings <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  pieces <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  collapse_unique_values(pieces)
}

max_or_na <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

lead_variant_from_pip <- function(variant_id, pip) {
  pip <- as.numeric(pip)
  if (length(pip) == 0 || all(is.na(pip))) return(NA_character_)
  collapse_unique_values(variant_id[pip == max(pip, na.rm = TRUE)])
}

## Standard display format: chromosome:position:allele1:allele2.
## This standardizes delimiters only; allele orientation is preserved from
## each validated input/membership key rather than re-inferred here.
## Handles individual IDs and semicolon-delimited fields.
standardize_variant <- function(x) {
  vapply(x, function(z) {
    if (length(z) == 0 || is.na(z) || z == "") return(NA_character_)
    
    vars <- unlist(strsplit(z, ";", fixed = TRUE), use.names = FALSE)
    vars <- trimws(vars)
    vars <- sub("^chr", "", vars)
    
    ## QTL CS representation: 11:118808595[b38]G,A
    vars <- sub(
      "^([^:]+):([0-9]+)\\[b38\\]([^,]+),(.+)$",
      "\\1:\\2:\\3:\\4",
      vars
    )
    
    ## ChromBPNet representation: chr10_324183_C_T or 10_324183_C_T
    vars <- sub(
      "^([^_]+)_([0-9]+)_([^_]+)_([^_]+)$",
      "\\1:\\2:\\3:\\4",
      vars
    )
    
    collapse_unique_values(vars)
  }, character(1))
}

combine_variant_fields <- function(x, y) {
  mapply(
    function(a, b) collapse_variant_strings(c(a, b)),
    x, y,
    USE.NAMES = FALSE
  )
}

count_variants_in_field <- function(x) {
  vapply(x, function(z) {
    if (length(z) == 0 || is.na(z) || z == "") return(0L)
    length(unique(unlist(strsplit(z, ";", fixed = TRUE), use.names = FALSE)))
  }, integer(1))
}

load_cs_files <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No CS files found in: ", dir, " using pattern: ", pattern)
  }
  
  message("Loading ", length(files), " CS files from: ", dir)
  print(basename(files))
  
  out <- data.table::rbindlist(
    lapply(files, data.table::fread),
    fill = TRUE
  )
  
  if ("cs" %in% names(out) && !"cs_index" %in% names(out)) {
    data.table::setnames(out, "cs", "cs_index")
  }
  
  out
}

## ======================================================================
## 2) LOAD CHROMBPNET VARIANTS AND DEFINE VALIDATED MEMBERSHIP SETS
## ======================================================================

cbpnet_vars <- fread(chrombpnet_file)

required_cbpnet_cols <- c("hg19_chr", "hg19_pos", "alt", "ref", "variant_id")
missing_cbpnet_cols <- setdiff(required_cbpnet_cols, names(cbpnet_vars))
if (length(missing_cbpnet_cols) > 0) {
  stop("ChromBPNet file is missing required columns: ", paste(missing_cbpnet_cols, collapse = ", "))
}

## Keep the validated GWAS/hg19 orientation from the working plotting workflow.
cbpnet_vars <- cbpnet_vars %>%
  mutate(
    hg19_chr = gsub("^chr", "", hg19_chr),
    hg19_var_lookup = paste0(hg19_chr, ":", hg19_pos, ":", alt, ":", ref)
  )

## Authoritative membership sets: do not replace these with reformatted lookups.
cbpnet_set37 <- unique(cbpnet_vars$hg19_var_lookup)
cbpnet_set_b38 <- unique(
  sub(
    pattern = "^chr([^_]+)_([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$",
    replacement = "\\1:\\2[b38]\\3,\\4",
    x = cbpnet_vars$variant_id
  )
)

## Display-only lookup from validated ChromBPNet variant fields.
cbpnet_lookup_display <- cbpnet_vars %>%
  transmute(
    chrombpnet_variant_hg37 = standardize_variant(hg19_var_lookup),
    chrombpnet_variant_hg38 = standardize_variant(variant_id),
    chrombpnet_variant_rsid = if ("rsid" %in% names(cbpnet_vars)) as.character(rsid) else NA_character_
  ) %>%
  distinct()

map_cbpnet_hg37_to_hg38 <- function(x) {
  vapply(x, function(z) {
    if (length(z) == 0 || is.na(z) || z == "") return(NA_character_)
    vars_hg37 <- unlist(strsplit(z, ";", fixed = TRUE), use.names = FALSE)
    cbpnet_lookup_display %>%
      filter(chrombpnet_variant_hg37 %in% vars_hg37) %>%
      pull(chrombpnet_variant_hg38) %>%
      collapse_unique_values()
  }, character(1))
}

## ======================================================================
## 3) LOAD MOLECULAR QTL CREDIBLE SETS AND ANNOTATE QTL-CS SUPPORT
## ======================================================================

## Load only the 22 chromosome-specific CS files. This avoids accidentally
## mixing in concatenated files and prevents partial-chromosome audit/parser
## mismatches.
eqtl_cs <- load_cs_files(
  eqtl_cs_dir,
  "^All_CD4T_cells_chr[0-9]+_credible_sets\\.txt$"
)

caqtl_cs <- load_cs_files(
  caqtl_cs_dir,
  "^CD4T_chromatin_chr[0-9]+_credible_sets\\.txt$"
)

if ("cs" %in% names(eqtl_cs)) data.table::setnames(eqtl_cs, "cs", "cs_index")
if ("cs" %in% names(caqtl_cs)) data.table::setnames(caqtl_cs, "cs", "cs_index")

required_eqtl_cols <- c("gene", "region", "cs_index", "variant_id")
required_caqtl_cols <- c("peak", "region", "cs_index", "variant_id")
if (length(setdiff(required_eqtl_cols, names(eqtl_cs))) > 0) {
  stop("eQTL CS data missing columns: ", paste(setdiff(required_eqtl_cols, names(eqtl_cs)), collapse = ", "))
}
if (length(setdiff(required_caqtl_cols, names(caqtl_cs))) > 0) {
  stop("caQTL CS data missing columns: ", paste(setdiff(required_caqtl_cols, names(caqtl_cs)), collapse = ", "))
}

## QC: confirm that the loaded CS universe spans chromosomes 1-22.
## This catches accidental partial loading, e.g. auditing against chr12-only CSs.
qtl_cs_chr_qc <- bind_rows(
  eqtl_cs %>%
    mutate(
      qtl_type = "eQTL",
      chr_from_variant = sub(":.*$", "", as.character(variant_id))
    ) %>%
    dplyr::count(qtl_type, chr_from_variant),
  
  caqtl_cs %>%
    mutate(
      qtl_type = "caQTL",
      chr_from_variant = sub(":.*$", "", as.character(variant_id))
    ) %>%
    dplyr::count(qtl_type, chr_from_variant)
)

message("QTL CS chromosome coverage:")
print(qtl_cs_chr_qc %>% arrange(qtl_type, as.integer(chr_from_variant)))

missing_qtl_cs_chr <- setdiff(as.character(1:22), unique(qtl_cs_chr_qc$chr_from_variant))
if (length(missing_qtl_cs_chr) > 0) {
  warning("Not all chromosomes 1-22 are represented in loaded QTL CS files: ",
          paste(missing_qtl_cs_chr, collapse = ","))
}

## This object is the authoritative QTL-CS ChromBPNet support call.
## It uses exactly the membership definition used in the original plotting workflow.
qtl_cs_cbpnet <- bind_rows(
  eqtl_cs %>%
    group_by(gene, region, cs_index) %>%
    summarise(
      has_cbpnet_qtl_cs = any(as.character(variant_id) %in% cbpnet_set_b38),
      n_qtl_variants = n(),
      .groups = "drop"
    ) %>%
    transmute(
      qtl_type = "eQTL",
      feature = as.character(gene),
      qtl_region = as.character(region),
      qtl_cs_index = as.integer(cs_index),
      has_cbpnet_qtl_cs,
      n_qtl_variants
    ),
  caqtl_cs %>%
    group_by(peak, region, cs_index) %>%
    summarise(
      has_cbpnet_qtl_cs = any(as.character(variant_id) %in% cbpnet_set_b38),
      n_qtl_variants = n(),
      .groups = "drop"
    ) %>%
    transmute(
      qtl_type = "caQTL",
      feature = as.character(peak),
      qtl_region = as.character(region),
      qtl_cs_index = as.integer(cs_index),
      has_cbpnet_qtl_cs,
      n_qtl_variants
    )
)

## ======================================================================
## 4) LOAD GWAS FINE-MAPPED CREDIBLE SETS
## ======================================================================

extract_gwas_cs_variants <- function(rds_path) {
  region_obj <- readRDS(rds_path)
  if (isFALSE(region_obj$converged)) return(NULL)
  
  cs_list <- region_obj$sets$cs
  if (is.null(cs_list) || length(cs_list) == 0) return(NULL)
  
  all_vars <- colnames(region_obj$alpha)
  if (is.null(all_vars)) stop("No variant column names found in alpha matrix: ", rds_path)
  
  bind_rows(lapply(seq_along(cs_list), function(i) {
    idx <- cs_list[[i]]
    data.frame(
      variant_id = all_vars[idx],
      cs_index = as.integer(i),
      pip = region_obj$pip[idx],
      stringsAsFactors = FALSE
    )
  }))
}

load_gwas_cs_variants <- function(gwas_id) {
  message("Loading GWAS CSs: ", gwas_id)
  
  files <- list.files(
    path = file.path(gwas_fm_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- files[grepl("_1\\.5Mb/susie/", files)]
  if (length(files) == 0) return(NULL)
  
  out <- bind_rows(lapply(files, function(f) {
    df <- extract_gwas_cs_variants(f)
    if (is.null(df)) return(NULL)
    df$region <- paste0("chr", sub(".*_([0-9]+\\.[0-9]+)\\.rds$", "\\1", f))
    df
  }))
  if (nrow(out) == 0) return(NULL)
  
  out %>%
    mutate(
      gwas_id = gwas_id,
      trait = extract_trait(gwas_id),
      region_cs = make_region_cs(region, cs_index),
      ## Authoritative GWAS-CS support call from original plotting workflow.
      has_cbpnet_gwas_variant = as.character(variant_id) %in% cbpnet_set37
    )
}

file_names <- unlist(
  fread(file.path(data.dir, "preprocessed_folder_names.txt"), header = FALSE),
  use.names = FALSE
)
file_names <- unique(as.character(file_names))

gwas_variant_all <- bind_rows(lapply(file_names, load_gwas_cs_variants))
if (nrow(gwas_variant_all) == 0) stop("No GWAS credible sets were loaded.")

gwas_cs_study <- gwas_variant_all %>%
  group_by(gwas_id, trait, region, cs_index, region_cs) %>%
  summarise(
    has_cbpnet_gwas_cs = any(has_cbpnet_gwas_variant),
    n_gwas_variants = n(),
    .groups = "drop"
  )

## ======================================================================
## 5) LOAD SIGNIFICANT COLOC PAIRS
## ======================================================================

coloc_eqtl_dir  <- file.path(data.dir, "coloc_results/eqtl_gwas_coloc")
coloc_caqtl_dir <- file.path(data.dir, "coloc_results")

resolve_qtl_region_col <- function(coloc_dt, qtl_cs, qtl_region_col = NULL, qtl_type) {
  if (!is.null(qtl_region_col)) {
    if (!qtl_region_col %in% names(coloc_dt)) {
      stop(qtl_type, " coloc file lacks column: ", qtl_region_col)
    }
    return(qtl_region_col)
  }
  
  candidates <- grep("region", names(coloc_dt), value = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0) stop(qtl_type, " coloc file has no region columns.")
  
  qtl_regions <- unique(as.character(qtl_cs$region))
  scores <- vapply(
    candidates,
    function(x) sum(unique(as.character(coloc_dt[[x]])) %in% qtl_regions),
    numeric(1)
  )
  
  if (all(scores == 0)) {
    stop(
      "Could not identify ", qtl_type, " QTL region column. Candidates: ",
      paste(candidates, collapse = ", ")
    )
  }
  
  chosen <- candidates[which.max(scores)]
  message("  Using ", qtl_type, " QTL region column: ", chosen)
  chosen
}

read_coloc_pairs <- function(path, gwas_id, qtl_type, feature_col, qtl_cs,
                             qtl_region_col = "region_GRC38") {
  if (!file.exists(path)) return(NULL)
  
  dt <- fread(path, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)
  
  needed <- c("pval", "pval_nominal", "region", "idx1", "idx2", feature_col, qtl_region_col)
  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) {
    stop(
      qtl_type, " coloc file missing: ", paste(missing, collapse = ", "),
      "\nFile: ", path,
      "\nAvailable columns: ", paste(names(dt), collapse = ", ")
    )
  }
  
  ## Coloc result files are assumed to have already been filtered at PP.H4.abf > 0.5
  ## during generation. Preserve PP.H4 when present, but do not re-filter here.
  pp4_candidates <- c("PP.H4.abf", "PP.H4", "PP.H4.abf.x", "PP.H4.abf.y")
  pp4_col <- pp4_candidates[pp4_candidates %in% names(dt)][1]
  
  if (length(pp4_col) == 0 || is.na(pp4_col)) {
    dt[, coloc_PP.H4.abf := NA_real_]
  } else {
    dt[, coloc_PP.H4.abf := as.numeric(get(pp4_col))]
    if (any(dt$coloc_PP.H4.abf <= PP_H4_THRESH, na.rm = TRUE)) {
      warning(
        qtl_type, " coloc file contains rows with PP.H4 <= ",
        PP_H4_THRESH, ": ", path,
        "\nFiles were expected to be pre-filtered."
      )
    }
  }
  
  ## Original association filters retained exactly.
  dt <- dt[pval < GWAS_P_THRESH & pval_nominal < QTL_P_THRESH]
  if (nrow(dt) == 0) return(NULL)
  
  ## Keep the coloc-reported QTL hit variant as descriptive metadata.
  ## ChromBPNet support is evaluated at the paired GWAS/QTL credible-set level,
  ## using idx1/idx2 to identify the colocated SuSiE components. The hit variant
  ## is not used to recover the QTL CS because it is selected from the coloc/LBF
  ## matrix, while the exported CS TSV contains only SuSiE CS members.
  qtl_hit_variant_candidates <- c("variant_id_GRC38", "eQTL_variant_GRC38", "caQTL_variant_GRC38")
  qtl_hit_variant_col <- qtl_hit_variant_candidates[qtl_hit_variant_candidates %in% names(dt)][1]
  
  if (length(qtl_hit_variant_col) == 0 || is.na(qtl_hit_variant_col)) {
    dt[, qtl_hit_variant_hg38 := NA_character_]
  } else {
    dt[, qtl_hit_variant_hg38 := as.character(get(qtl_hit_variant_col))]
  }
  
  ## IMPORTANT:
  ## For GWAS-QTL coloc outputs, idx1 and idx2 are interpreted as the
  ## GWAS and QTL SuSiE CS/component indices, respectively. The paired
  ## QTL CS is recovered using:
  ##   qtl_type + feature + region_GRC38 + idx2
  ## This reports confident CS-level colocalization and ChromBPNet membership.
  ## The coloc hit variant is reported separately and should not be treated as
  ## the sole fine-mapped QTL variant unless it is also present in the CS TSV.
  out <- dt %>%
    as_tibble() %>%
    transmute(
      gwas_id = gwas_id,
      trait = extract_trait(gwas_id),
      
      ## GWAS CS side
      region = as.character(region),
      cs_index = as.integer(idx1),
      region_cs = make_region_cs(region, idx1),
      
      ## QTL CS side
      qtl_type = qtl_type,
      feature = as.character(.data[[feature_col]]),
      qtl_region = as.character(.data[[qtl_region_col]]),
      qtl_cs_index = as.integer(idx2),
      
      ## Descriptive coloc metadata, not used for CS recovery
      qtl_hit_variant_hg38 = as.character(qtl_hit_variant_hg38),
      coloc_PP.H4.abf = coloc_PP.H4.abf,
      coloc_source_file = path
    ) %>%
    distinct()
  
  ## Per-file QC: parsed QTL CS keys should exist in the corresponding
  ## fine-mapping CS table.
  qtl_feature_lookup_col <- if (qtl_type == "eQTL") "gene" else "peak"
  
  qtl_key_check <- out %>%
    distinct(qtl_type, feature, qtl_region, qtl_cs_index) %>%
    left_join(
      qtl_cs %>%
        transmute(
          feature = as.character(.data[[qtl_feature_lookup_col]]),
          qtl_region = as.character(region),
          qtl_cs_index = as.integer(cs_index),
          qtl_cs_lookup_matched = TRUE
        ) %>%
        distinct(),
      by = c("feature", "qtl_region", "qtl_cs_index")
    ) %>%
    mutate(qtl_cs_lookup_matched = replace_na(qtl_cs_lookup_matched, FALSE))
  
  if (any(!qtl_key_check$qtl_cs_lookup_matched)) {
    warning(
      qtl_type, " coloc rows parsed from ", basename(path),
      " include QTL CS keys not found in the loaded fine-mapping CS table. ",
      "This usually indicates wrong CS files or region/feature naming mismatch."
    )
  }
  
  out
}

load_study_coloc_pairs <- function(gwas_id) {
  eqtl_path <- file.path(coloc_eqtl_dir, paste0(gwas_id, "_coloc_results.csv"))
  caqtl_path <- file.path(
    coloc_caqtl_dir,
    paste0("CD4T_chromatin_", gwas_id, "_coloc_results.csv")
  )
  
  bind_rows(
    read_coloc_pairs(
      eqtl_path, gwas_id, "eQTL", "gene", eqtl_cs,
      qtl_region_col = "region_GRC38"
    ),
    read_coloc_pairs(
      caqtl_path, gwas_id, "caQTL", "peak", caqtl_cs,
      qtl_region_col = "region_GRC38"
    )
  )
}

coloc_pairs_raw <- bind_rows(lapply(file_names, load_study_coloc_pairs)) %>%
  distinct()

## Global QC: parsed QTL-CS keys should match the QTL-CS support lookup.
## This checks the CS-based parser: qtl_type + feature + region_GRC38 + idx2.
coloc_qtl_key_qc <- coloc_pairs_raw %>%
  distinct(qtl_type, feature, qtl_region, qtl_cs_index) %>%
  left_join(
    qtl_cs_cbpnet %>%
      select(all_of(KEY_QTL), n_qtl_variants) %>%
      distinct(),
    by = KEY_QTL
  ) %>%
  mutate(qtl_cs_key_matched = !is.na(n_qtl_variants))

message("GWAS-QTL coloc parsed QTL-CS key match summary:")
print(
  coloc_qtl_key_qc %>%
    dplyr::count(qtl_type, qtl_cs_key_matched)
)

unmatched_coloc_qtl_keys <- coloc_qtl_key_qc %>%
  filter(!qtl_cs_key_matched)

write.csv(
  unmatched_coloc_qtl_keys,
  file.path(output_dir, "QC_unmatched_coloc_qtl_keys_after_parsing.csv"),
  row.names = FALSE
)

if (nrow(unmatched_coloc_qtl_keys) > 0) {
  warning(
    nrow(unmatched_coloc_qtl_keys),
    " parsed QTL-CS keys from coloc results do not match qtl_cs_cbpnet. ",
    "Inspect QC_unmatched_coloc_qtl_keys_after_parsing.csv."
  )
}

## ======================================================================
## 6) PRIMARY STUDY-SPECIFIC COLOC UNIVERSE
## ======================================================================
##
## A raw significant coloc pair enters the analysis only when its matching
## GWAS credible set exists in the same GWAS study. This is the explicit
## equivalent of the original plotting logic, which joined study-specific
## coloc annotations onto study-specific GWAS CSs before deduplicating loci.

## `distinct(across())` can fail in older dplyr versions; select + distinct is stable.
gwas_study_key <- gwas_cs_study %>%
  select(all_of(KEY_STUDY)) %>%
  distinct()

excluded_coloc_pairs_not_in_primary_gwas_cs <- coloc_pairs_raw %>%
  anti_join(gwas_study_key, by = KEY_STUDY) %>%
  arrange(trait, region, cs_index, gwas_id, qtl_type, feature)

message(
  "Significant coloc feature rows excluded because their study-specific GWAS CS ",
  "was absent from the primary GWAS-CS universe: ",
  nrow(excluded_coloc_pairs_not_in_primary_gwas_cs)
)

coloc_pairs_primary_study <- coloc_pairs_raw %>%
  inner_join(gwas_study_key, by = KEY_STUDY) %>%
  left_join(qtl_cs_cbpnet, by = KEY_QTL) %>%
  mutate(
    qtl_cs_matched = !is.na(n_qtl_variants),
    ## Same behavior as the plot workflow: an unmatched QTL CS does not
    ## provide ChromBPNet support, but does not erase GWAS-CS support.
    has_cbpnet_qtl_cs = replace_na(has_cbpnet_qtl_cs, FALSE)
  )

unmatched_primary_pairs_to_qtl_cs <- coloc_pairs_primary_study %>%
  filter(!qtl_cs_matched) %>%
  arrange(trait, region, cs_index, gwas_id, qtl_type, feature)

message(
  "Primary-analysis coloc feature rows without a matched QTL CS lookup: ",
  nrow(unmatched_primary_pairs_to_qtl_cs)
)

## ======================================================================
## 7) CANONICAL LOCUS-LEVEL OBJECT: USED BY BOTH PLOTS AND TABLE
## ======================================================================

coloc_cs_study <- coloc_pairs_primary_study %>%
  group_by(gwas_id, trait, region, cs_index, region_cs) %>%
  summarise(
    eqtl_coloc = any(qtl_type == "eQTL"),
    caqtl_coloc = any(qtl_type == "caQTL"),
    has_cbpnet_qtl_cs = any(has_cbpnet_qtl_cs),
    n_paired_qtl_cs = n_distinct(
      paste(qtl_type, feature, qtl_region, qtl_cs_index, sep = "|")
    ),
    .groups = "drop"
  )

gwas_cs_study_annot <- gwas_cs_study %>%
  left_join(coloc_cs_study, by = KEY_STUDY) %>%
  mutate(
    eqtl_coloc = replace_na(eqtl_coloc, FALSE),
    caqtl_coloc = replace_na(caqtl_coloc, FALSE),
    has_cbpnet_qtl_cs = replace_na(has_cbpnet_qtl_cs, FALSE),
    n_paired_qtl_cs = replace_na(n_paired_qtl_cs, 0L)
  )

trait_locus_master <- gwas_cs_study_annot %>%
  group_by(trait, region, cs_index, region_cs) %>%
  summarise(
    n_studies_at_locus = n_distinct(gwas_id),
    study_ids = paste(sort(unique(gwas_id)), collapse = ";"),
    has_cbpnet_gwas_cs = any(has_cbpnet_gwas_cs),
    eqtl_coloc = any(eqtl_coloc),
    caqtl_coloc = any(caqtl_coloc),
    has_cbpnet_qtl_cs = any(has_cbpnet_qtl_cs),
    n_paired_qtl_cs_across_studies = sum(n_paired_qtl_cs),
    .groups = "drop"
  ) %>%
  mutate(
    qtl_coloc_cs = eqtl_coloc | caqtl_coloc,
    
    qtl_coloc_category = case_when(
      eqtl_coloc & caqtl_coloc  ~ "both",
      eqtl_coloc & !caqtl_coloc ~ "just_eqtl",
      !eqtl_coloc & caqtl_coloc ~ "just_caqtl",
      TRUE                      ~ "no_coloc"
    ),
    
    has_cbpnet_either_cs = has_cbpnet_gwas_cs | has_cbpnet_qtl_cs,
    
    cbpnet_membership = case_when(
      qtl_coloc_cs & has_cbpnet_gwas_cs & has_cbpnet_qtl_cs ~ "Both_GWAS_and_QTL_CS",
      qtl_coloc_cs & has_cbpnet_gwas_cs & !has_cbpnet_qtl_cs ~ "GWAS_CS_only_coloc",
      qtl_coloc_cs & !has_cbpnet_gwas_cs & has_cbpnet_qtl_cs ~ "QTL_CS_only_coloc",
      !qtl_coloc_cs & has_cbpnet_gwas_cs ~ "GWAS_CS_only_noncoloc",
      TRUE ~ "No_ChromBPNet"
    ),
    
    cbpnet_category = case_when(
      qtl_coloc_cs & has_cbpnet_either_cs ~ "QTL_coloc_with_ChromBPNet",
      !qtl_coloc_cs & has_cbpnet_gwas_cs ~ "ChromBPNet_specific",
      qtl_coloc_cs & !has_cbpnet_either_cs ~ "QTL_coloc_no_ChromBPNet",
      TRUE ~ "Neither"
    )
  )

stopifnot(
  nrow(
    trait_locus_master %>%
      count(trait, region_cs) %>%
      filter(n > 1)
  ) == 0
)

repeated_trait_loci <- trait_locus_master %>%
  filter(n_studies_at_locus > 1) %>%
  arrange(desc(n_studies_at_locus), trait, region_cs)

message("Repeated trait-loci collapsed across studies: ", nrow(repeated_trait_loci))

## Back-compatible alias, if other downstream code uses the former name.
trait_locus_summary <- trait_locus_master

## ======================================================================
## 8) PLOT SUMMARIES FROM THE CANONICAL LOCUS TABLE
## ======================================================================

coloc_table_summary <- trait_locus_master %>%
  group_by(trait) %>%
  summarise(
    total_cond_indep = n(),
    just_eqtl = sum(qtl_coloc_category == "just_eqtl"),
    just_caqtl = sum(qtl_coloc_category == "just_caqtl"),
    both = sum(qtl_coloc_category == "both"),
    n_repeated_loci_collapsed = sum(n_studies_at_locus > 1),
    .groups = "drop"
  ) %>%
  mutate(
    trait_label = label_trait(trait),
    total_coloc = just_eqtl + just_caqtl + both,
    frac_coloc = ifelse(total_cond_indep > 0, total_coloc / total_cond_indep, 0)
  )

trait_summary <- trait_locus_master %>%
  group_by(trait) %>%
  summarise(
    total_cs_all = n(),
    cs_with_cbpnet_in_gwas = sum(has_cbpnet_gwas_cs),
    cs_with_cbpnet_support = sum(has_cbpnet_either_cs),
    cs_cbpnet_specific = sum(cbpnet_category == "ChromBPNet_specific"),
    cs_qtl_coloc_with_cbpnet = sum(cbpnet_category == "QTL_coloc_with_ChromBPNet"),
    cs_qtl_coloc_total = sum(qtl_coloc_cs),
    cs_qtl_coloc_no_cbpnet = sum(cbpnet_category == "QTL_coloc_no_ChromBPNet"),
    cs_coloc_gwas_only = sum(cbpnet_membership == "GWAS_CS_only_coloc"),
    cs_coloc_qtl_only = sum(cbpnet_membership == "QTL_CS_only_coloc"),
    cs_coloc_both_sets = sum(cbpnet_membership == "Both_GWAS_and_QTL_CS"),
    n_repeated_loci_collapsed = sum(n_studies_at_locus > 1),
    .groups = "drop"
  ) %>%
  mutate(
    trait_label = label_trait(trait),
    prop_cs_with_cbpnet_support = ifelse(total_cs_all > 0, cs_with_cbpnet_support / total_cs_all, 0),
    prop_cbpnet_in_gwas = ifelse(total_cs_all > 0, cs_with_cbpnet_in_gwas / total_cs_all, 0),
    prop_cbpnet_specific = ifelse(total_cs_all > 0, cs_cbpnet_specific / total_cs_all, 0),
    prop_qtl_coloc_with_cbpnet = ifelse(total_cs_all > 0, cs_qtl_coloc_with_cbpnet / total_cs_all, 0)
  )

denominator_qc <- coloc_table_summary %>%
  select(trait, total_plot2 = total_cond_indep) %>%
  left_join(trait_summary %>% select(trait, total_plot3 = total_cs_all), by = "trait") %>%
  mutate(passes_check = total_plot2 == total_plot3)

category_qc <- trait_summary %>%
  mutate(
    category_sum = cs_cbpnet_specific + cs_qtl_coloc_with_cbpnet,
    passes_check = cs_with_cbpnet_support == category_sum
  )

stopifnot(all(denominator_qc$passes_check), all(category_qc$passes_check))

## Only display traits that have at least one QTL-GWAS colocalized locus.
keep_traits <- coloc_table_summary %>%
  filter(total_coloc > 0) %>%
  pull(trait)

coloc_table_summary_plot <- coloc_table_summary %>% filter(trait %in% keep_traits)
trait_summary_plot <- trait_summary %>% filter(trait %in% keep_traits)

trait_order <- coloc_table_summary_plot %>%
  arrange(total_cond_indep, trait_label) %>%
  pull(trait_label)

## ======================================================================
## 9) PLOT 2: QTL-GWAS COLOCALIZATION AT UNIQUE TRAIT LOCI
## ======================================================================

plot_df <- coloc_table_summary_plot %>%
  select(trait, trait_label, total_cond_indep, just_eqtl, just_caqtl, both) %>%
  pivot_longer(
    cols = c(just_eqtl, just_caqtl, both),
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(
    proportion = ifelse(total_cond_indep > 0, count / total_cond_indep, 0),
    trait_label = factor(trait_label, levels = trait_order)
  )

avg_total_coloc <- coloc_table_summary_plot %>%
  select(trait_label, total_cond_indep) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order)) %>%
  arrange(total_cond_indep)

max_total_coloc <- max(avg_total_coloc$total_cond_indep, na.rm = TRUE)
scaling_factor_coloc <- ifelse(max_total_coloc > 0, 1 / max_total_coloc, 1)

plot2 <- ggplot(plot_df, aes(x = trait_label, y = proportion, fill = category)) +
  geom_col(color = "gray40", linewidth = 0.1) +
  geom_text(
    data = subset(plot_df, count != 0),
    aes(label = count),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6
  ) +
  geom_line(
    data = avg_total_coloc,
    aes(x = trait_label, y = total_cond_indep * scaling_factor_coloc, group = 1),
    inherit.aes = FALSE,
    color = "gray50",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    data = avg_total_coloc,
    aes(x = trait_label, y = total_cond_indep * scaling_factor_coloc),
    inherit.aes = FALSE,
    shape = 21,
    fill = "gray30",
    size = 4,
    stroke = 0.3,
    alpha = 0.5
  ) +
  scale_fill_manual(values = qtl_cmap_coloc) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0.01, 0),
    labels = percent_format(accuracy = 1),
    name = "Proportion of unique GWAS loci colocalized with QTLs",
    sec.axis = sec_axis(
      ~ . / scaling_factor_coloc,
      breaks = pretty(c(0, max_total_coloc)),
      name = "Total unique GWAS loci per trait"
    )
  ) +
  coord_flip() +
  labs(title = "Summary of QTL-GWAS colocalization", x = "GWAS trait") +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x.bottom = element_text(size = 12, color = "black"),
    axis.text.x.top = element_text(size = 12, color = "gray50"),
    axis.title.x.top = element_text(size = 14, color = "gray50")
  )

## ======================================================================
## 10) PLOT 3: CHROMBPNET SUPPORT AT UNIQUE TRAIT LOCI
## ======================================================================

trait_long <- trait_summary_plot %>%
  select(trait, trait_label, total_cs_all, cs_cbpnet_specific, cs_qtl_coloc_with_cbpnet) %>%
  pivot_longer(
    cols = c(cs_cbpnet_specific, cs_qtl_coloc_with_cbpnet),
    names_to = "cbpnet_cat",
    values_to = "cs_n"
  ) %>%
  mutate(
    cbpnet_cat = recode(
      cbpnet_cat,
      cs_cbpnet_specific = "ChromBPNet_specific",
      cs_qtl_coloc_with_cbpnet = "QTL_coloc_with_ChromBPNet"
    ),
    cbpnet_cat = factor(cbpnet_cat, levels = c("ChromBPNet_specific", "QTL_coloc_with_ChromBPNet")),
    prop = ifelse(total_cs_all > 0, cs_n / total_cs_all, 0),
    trait_label = factor(trait_label, levels = trait_order)
  )

trait_summary_lab <- trait_summary_plot %>%
  filter(cs_with_cbpnet_support > 0) %>%
  mutate(trait_label = factor(trait_label, levels = trait_order))

plot3 <- ggplot(trait_long, aes(x = trait_label, y = prop, fill = cbpnet_cat)) +
  geom_col(color = "gray40", linewidth = 0.1) +
  geom_text(
    data = trait_summary_lab,
    aes(x = trait_label, y = prop_cs_with_cbpnet_support / 2, label = cs_with_cbpnet_support),
    inherit.aes = FALSE,
    color = "black",
    size = 6
  ) +
  scale_fill_manual(
    values = c(
      "ChromBPNet_specific" = "#0071ba",
      "QTL_coloc_with_ChromBPNet" = "#f7d58b"
    ),
    labels = c(
      "ChromBPNet_specific" = "ChromBPNet only",
      "QTL_coloc_with_ChromBPNet" = "ChromBPNet + QTL colocalization"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 0.30),
    breaks = seq(0, 0.30, 0.05),
    labels = percent_format(accuracy = 1),
    expand = c(0.01, 0),
    name = "Proportion of unique GWAS loci with ChromBPNet support"
  ) +
  coord_flip() +
  theme_cowplot() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

plot2_compact <- plot2 + theme(plot.margin = margin(t = 5, r = 2, b = 5, l = 5))
plot3_compact <- plot3 + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 2))

combined_plot <- cowplot::plot_grid(
  plot2_compact,
  plot3_compact,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(2, 0.5)
)

ggsave(
  filename = plot_output_file,
  plot = combined_plot,
  height = 17,
  width = 17,
  dpi = 300,
  units = "in"
)

## ======================================================================
## 11) DESCRIPTIVE GWAS-CS DETAIL FOR SUPPLEMENTARY TABLE 8
## ======================================================================
##
## These columns describe retained loci; they do not classify locus inclusion.

gwas_cs_detail_study <- gwas_variant_all %>%
  mutate(
    variant_id_hg37 = standardize_variant(variant_id)
  ) %>%
  group_by(gwas_id, trait, region, cs_index, region_cs) %>%
  summarise(
    gwas_lead_variant_hg37 = lead_variant_from_pip(variant_id_hg37, pip),
    gwas_lead_pip = max_or_na(pip),
    gwas_cs_variants_hg37 = collapse_unique_values(variant_id_hg37),
    chrombpnet_variants_gwas_cs_hg37 = collapse_unique_values(
      variant_id_hg37[has_cbpnet_gwas_variant]
    ),
    detail_has_cbpnet_gwas_cs = any(has_cbpnet_gwas_variant),
    .groups = "drop"
  )

gwas_cs_detail_trait_locus <- gwas_cs_detail_study %>%
  group_by(trait, region, cs_index, region_cs) %>%
  summarise(
    n_gwas_studies_at_locus = n_distinct(gwas_id),
    contributing_gwas_studies = collapse_unique_values(gwas_id),
    gwas_lead_variants_hg37 = collapse_unique_values(gwas_lead_variant_hg37),
    gwas_lead_variants_by_study_hg37 = collapse_unique_values(
      paste0(gwas_id, "=", gwas_lead_variant_hg37)
    ),
    maximum_gwas_lead_pip = max_or_na(gwas_lead_pip),
    gwas_cs_variants_hg37 = collapse_variant_strings(gwas_cs_variants_hg37),
    chrombpnet_variants_gwas_cs_hg37 = collapse_variant_strings(
      chrombpnet_variants_gwas_cs_hg37
    ),
    detail_has_cbpnet_gwas_cs = any(detail_has_cbpnet_gwas_cs),
    .groups = "drop"
  ) %>%
  mutate(
    chrombpnet_variants_gwas_cs_hg38 = map_cbpnet_hg37_to_hg38(
      chrombpnet_variants_gwas_cs_hg37
    ),
    n_chrombpnet_variants_gwas_cs = count_variants_in_field(
      chrombpnet_variants_gwas_cs_hg37
    )
  )

gwas_support_qc <- trait_locus_master %>%
  select(all_of(KEY_LOCUS), master_has_cbpnet_gwas_cs = has_cbpnet_gwas_cs) %>%
  left_join(
    gwas_cs_detail_trait_locus %>%
      select(all_of(KEY_LOCUS), detail_has_cbpnet_gwas_cs),
    by = KEY_LOCUS
  ) %>%
  mutate(passes_check = master_has_cbpnet_gwas_cs == detail_has_cbpnet_gwas_cs)

stopifnot(all(gwas_support_qc$passes_check))

## ======================================================================
## 12) DESCRIPTIVE QTL-CS DETAIL USING THE SAME SUPPORT CALL AS THE MASTER
## ======================================================================

summarise_qtl_cs_detail <- function(qtl_dt, qtl_type_value, feature_col) {
  pip_candidates <- c("pip", "PIP", "posterior_inclusion_probability", "PIP.y")
  pip_col <- pip_candidates[pip_candidates %in% names(qtl_dt)][1]
  
  dat <- as_tibble(qtl_dt) %>%
    mutate(
      feature = as.character(.data[[feature_col]]),
      qtl_region = as.character(region),
      qtl_cs_index = as.integer(cs_index),
      original_variant_id_hg38 = as.character(variant_id),
      variant_id_hg38 = standardize_variant(original_variant_id_hg38),
      ## Critical: exact Boolean support definition used above and in the plots.
      is_cbpnet_qtl_variant = original_variant_id_hg38 %in% cbpnet_set_b38
    )
  
  if (length(pip_col) == 0 || is.na(pip_col)) {
    warning(
      "No PIP column detected in ", qtl_type_value,
      " CS table; qtl_lead_variant_hg38 and qtl_lead_pip are NA."
    )
    
    out <- dat %>%
      group_by(feature, qtl_region, qtl_cs_index) %>%
      summarise(
        qtl_lead_variant_hg38 = NA_character_,
        qtl_lead_pip = NA_real_,
        qtl_cs_variants_hg38 = collapse_unique_values(variant_id_hg38),
        chrombpnet_variants_qtl_cs_hg38 = collapse_unique_values(
          variant_id_hg38[is_cbpnet_qtl_variant]
        ),
        detail_has_cbpnet_qtl_cs = any(is_cbpnet_qtl_variant),
        .groups = "drop"
      )
  } else {
    out <- dat %>%
      group_by(feature, qtl_region, qtl_cs_index) %>%
      summarise(
        qtl_lead_variant_hg38 = lead_variant_from_pip(
          variant_id_hg38, .data[[pip_col]]
        ),
        qtl_lead_pip = max_or_na(.data[[pip_col]]),
        qtl_cs_variants_hg38 = collapse_unique_values(variant_id_hg38),
        chrombpnet_variants_qtl_cs_hg38 = collapse_unique_values(
          variant_id_hg38[is_cbpnet_qtl_variant]
        ),
        detail_has_cbpnet_qtl_cs = any(is_cbpnet_qtl_variant),
        .groups = "drop"
      )
  }
  
  out %>%
    mutate(
      qtl_type = qtl_type_value,
      n_chrombpnet_variants_qtl_cs = count_variants_in_field(
        chrombpnet_variants_qtl_cs_hg38
      )
    ) %>%
    select(
      qtl_type,
      feature,
      qtl_region,
      qtl_cs_index,
      qtl_lead_variant_hg38,
      qtl_lead_pip,
      qtl_cs_variants_hg38,
      chrombpnet_variants_qtl_cs_hg38,
      n_chrombpnet_variants_qtl_cs,
      detail_has_cbpnet_qtl_cs
    )
}

qtl_cs_detail <- bind_rows(
  summarise_qtl_cs_detail(eqtl_cs, "eQTL", "gene"),
  summarise_qtl_cs_detail(caqtl_cs, "caQTL", "peak")
)

qtl_support_qc <- qtl_cs_cbpnet %>%
  select(all_of(KEY_QTL), master_has_cbpnet_qtl_cs = has_cbpnet_qtl_cs) %>%
  full_join(
    qtl_cs_detail %>%
      select(all_of(KEY_QTL), detail_has_cbpnet_qtl_cs),
    by = KEY_QTL
  ) %>%
  mutate(
    master_has_cbpnet_qtl_cs = replace_na(master_has_cbpnet_qtl_cs, FALSE),
    detail_has_cbpnet_qtl_cs = replace_na(detail_has_cbpnet_qtl_cs, FALSE),
    passes_check = master_has_cbpnet_qtl_cs == detail_has_cbpnet_qtl_cs
  )

if (any(!qtl_support_qc$passes_check)) {
  write.csv(
    qtl_support_qc %>% filter(!passes_check),
    file.path(output_dir, "QC_qtl_chrombpnet_support_disagreement.csv"),
    row.names = FALSE
  )
  stop(
    "QTL ChromBPNet support calls in descriptive details do not match the ",
    "canonical plot/table support call. Inspect QC_qtl_chrombpnet_support_disagreement.csv."
  )
}

## ======================================================================
## 13) LOCUS-LEVEL QTL PAIRING DETAILS FOR SUPPLEMENTARY TABLE 8
## ======================================================================
##
## The plot is defined at the unique trait-GWAS-locus level:
##   trait + region + cs_index + region_cs
##
## Supplementary Table 8 will now use the same unit.
## Multiple colocalized QTL credible sets at the same locus are collapsed
## into semicolon-separated descriptive fields rather than separate rows.

## ----------------------------------------------------------------------
## 13A) Additional helper functions
## ----------------------------------------------------------------------

## Collapse fields that may already contain semicolon-separated values.
collapse_semicolon_fields <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  pieces <- unlist(
    strsplit(x, ";", fixed = TRUE),
    use.names = FALSE
  )
  
  pieces <- trimws(pieces)
  collapse_unique_values(pieces)
}


## Count unique values across fields that may already contain
## semicolon-separated entries.
count_semicolon_values <- function(x) {
  collapsed <- collapse_semicolon_fields(x)
  
  if (is.na(collapsed) || collapsed == "") {
    return(0L)
  }
  
  as.integer(
    length(
      unique(
        unlist(
          strsplit(collapsed, ";", fixed = TRUE),
          use.names = FALSE
        )
      )
    )
  )
}


## Construct and collapse QTL credible-set identifiers safely.
collapse_qtl_cs_keys <- function(
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index
) {
  keep <- !is.na(qtl_type) &
    qtl_type != "" &
    !is.na(feature) &
    feature != "" &
    !is.na(qtl_region) &
    qtl_region != "" &
    !is.na(qtl_cs_index)
  
  if (!any(keep)) {
    return(NA_character_)
  }
  
  keys <- paste(
    qtl_type[keep],
    feature[keep],
    qtl_region[keep],
    qtl_cs_index[keep],
    sep = ":"
  )
  
  collapse_unique_values(keys)
}


## Count unique QTL credible-set identifiers safely.
count_qtl_cs_keys <- function(
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index
) {
  keep <- !is.na(qtl_type) &
    qtl_type != "" &
    !is.na(feature) &
    feature != "" &
    !is.na(qtl_region) &
    qtl_region != "" &
    !is.na(qtl_cs_index)
  
  if (!any(keep)) {
    return(0L)
  }
  
  as.integer(
    n_distinct(
      paste(
        qtl_type[keep],
        feature[keep],
        qtl_region[keep],
        qtl_cs_index[keep],
        sep = "|"
      )
    )
  )
}


## ----------------------------------------------------------------------
## 13B) Create one row per unique GWAS-locus–QTL-CS pairing
## ----------------------------------------------------------------------

coloc_pair_trait_locus <- coloc_pairs_primary_study %>%
  select(
    gwas_id,
    trait,
    region,
    cs_index,
    region_cs,
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index,
    qtl_hit_variant_hg38,
    coloc_PP.H4.abf,
    has_cbpnet_qtl_cs
  ) %>%
  distinct() %>%
  group_by(
    trait,
    region,
    cs_index,
    region_cs,
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index
  ) %>%
  summarise(
    n_coloc_supporting_studies = n_distinct(gwas_id),
    
    coloc_supporting_gwas_studies =
      collapse_unique_values(gwas_id),
    
    coloc_qtl_hit_variants_hg38 =
      collapse_unique_values(qtl_hit_variant_hg38),
    
    max_coloc_PP.H4.abf =
      max_or_na(coloc_PP.H4.abf),
    
    paired_qtl_cs_has_cbpnet =
      any(has_cbpnet_qtl_cs),
    
    .groups = "drop"
  )

## ----------------------------------------------------------------------
## Define the exact supported QTL-colocalized locus universe from Plot 3
## ----------------------------------------------------------------------

supported_coloc_loci <- trait_locus_master %>%
  filter(
    qtl_coloc_cs,
    has_cbpnet_either_cs
  ) %>%
  transmute(
    trait,
    region,
    cs_index,
    region_cs,
    
    qtl_coloc_category,
    
    locus_has_cbpnet_gwas_cs =
      has_cbpnet_gwas_cs,
    
    locus_has_cbpnet_qtl_cs =
      has_cbpnet_qtl_cs,
    
    locus_has_cbpnet_either_cs =
      has_cbpnet_either_cs
  )

## One row per unique trait-GWAS locus
supported_coloc_locus_duplicates <- supported_coloc_loci %>%
  count(
    trait,
    region,
    cs_index,
    region_cs,
    name = "n"
  ) %>%
  filter(n > 1)

if (nrow(supported_coloc_locus_duplicates) > 0) {
  stop(
    "supported_coloc_loci contains duplicated trait-GWAS loci."
  )
}

coloc_pair_detail_trait_locus <- coloc_pair_trait_locus %>%
  left_join(
    qtl_cs_detail,
    by = KEY_QTL
  ) %>%
  mutate(
    qtl_cs_detail_matched = !is.na(qtl_cs_variants_hg38),
    
    ## Use the authoritative pairing-level support call.
    QTL_CS_has_ChromBPNet_variant = paired_qtl_cs_has_cbpnet,
    
    qtl_cs_key = paste(
      qtl_type,
      feature,
      qtl_region,
      qtl_cs_index,
      sep = "|"
    )
  )


## QC: each locus-QTL-CS combination should occur only once.
duplicate_locus_qtl_pairs <- coloc_pair_detail_trait_locus %>%
  count(
    across(all_of(KEY_LOCUS)),
    qtl_type,
    feature,
    qtl_region,
    qtl_cs_index,
    name = "n"
  ) %>%
  filter(n > 1)

if (nrow(duplicate_locus_qtl_pairs) > 0) {
  write.csv(
    duplicate_locus_qtl_pairs,
    file.path(
      output_dir,
      "QC_duplicate_locus_QTL_CS_pairings.csv"
    ),
    row.names = FALSE
  )
  
  stop(
    "Duplicate locus-QTL credible-set pairings remain before locus-level ",
    "collapse. Inspect QC_duplicate_locus_QTL_CS_pairings.csv."
  )
}


## ----------------------------------------------------------------------
## 13C) Collapse all paired QTL information to one row per GWAS locus
## ----------------------------------------------------------------------

qtl_pair_summary_at_locus <- coloc_pair_detail_trait_locus %>%
  group_by(across(all_of(KEY_LOCUS))) %>%
  summarise(
    ## ---------------------------------------------------------------
    ## All colocalized molecular features
    ## ---------------------------------------------------------------
    
    linked_eGenes = collapse_unique_values(
      feature[qtl_type == "eQTL"]
    ),
    
    linked_caPeaks = collapse_unique_values(
      feature[qtl_type == "caQTL"]
    ),
    
    n_linked_eGenes = n_distinct(
      feature[
        qtl_type == "eQTL" &
          !is.na(feature) &
          feature != ""
      ]
    ),
    
    n_linked_caPeaks = n_distinct(
      feature[
        qtl_type == "caQTL" &
          !is.na(feature) &
          feature != ""
      ]
    ),
    
    ## ---------------------------------------------------------------
    ## All paired QTL credible sets
    ## ---------------------------------------------------------------
    
    all_paired_qtl_credible_sets = collapse_qtl_cs_keys(
      qtl_type,
      feature,
      qtl_region,
      qtl_cs_index
    ),
    
    n_paired_qtl_credible_sets = count_qtl_cs_keys(
      qtl_type,
      feature,
      qtl_region,
      qtl_cs_index
    ),
    
    paired_eQTL_credible_sets = collapse_qtl_cs_keys(
      qtl_type[qtl_type == "eQTL"],
      feature[qtl_type == "eQTL"],
      qtl_region[qtl_type == "eQTL"],
      qtl_cs_index[qtl_type == "eQTL"]
    ),
    
    paired_caQTL_credible_sets = collapse_qtl_cs_keys(
      qtl_type[qtl_type == "caQTL"],
      feature[qtl_type == "caQTL"],
      qtl_region[qtl_type == "caQTL"],
      qtl_cs_index[qtl_type == "caQTL"]
    ),
    
    ## ---------------------------------------------------------------
    ## Only paired QTL CSs that themselves contain ChromBPNet support
    ## ---------------------------------------------------------------
    
    chrombpnet_supporting_qtl_credible_sets = collapse_qtl_cs_keys(
      qtl_type[QTL_CS_has_ChromBPNet_variant],
      feature[QTL_CS_has_ChromBPNet_variant],
      qtl_region[QTL_CS_has_ChromBPNet_variant],
      qtl_cs_index[QTL_CS_has_ChromBPNet_variant]
    ),
    
    n_chrombpnet_supporting_qtl_credible_sets = count_qtl_cs_keys(
      qtl_type[QTL_CS_has_ChromBPNet_variant],
      feature[QTL_CS_has_ChromBPNet_variant],
      qtl_region[QTL_CS_has_ChromBPNet_variant],
      qtl_cs_index[QTL_CS_has_ChromBPNet_variant]
    ),
    
    chrombpnet_supporting_eQTL_credible_sets = collapse_qtl_cs_keys(
      qtl_type[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "eQTL"
      ],
      feature[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "eQTL"
      ],
      qtl_region[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "eQTL"
      ],
      qtl_cs_index[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "eQTL"
      ]
    ),
    
    chrombpnet_supporting_caQTL_credible_sets = collapse_qtl_cs_keys(
      qtl_type[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "caQTL"
      ],
      feature[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "caQTL"
      ],
      qtl_region[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "caQTL"
      ],
      qtl_cs_index[
        QTL_CS_has_ChromBPNet_variant &
          qtl_type == "caQTL"
      ]
    ),
    
    ## ---------------------------------------------------------------
    ## QTL-CS fine-mapping details
    ## ---------------------------------------------------------------
    
    qtl_lead_variants_hg38 = collapse_variant_strings(
      qtl_lead_variant_hg38
    ),
    
    maximum_qtl_lead_pip = max_or_na(
      qtl_lead_pip
    ),
    
    all_paired_qtl_cs_variants_hg38 = collapse_variant_strings(
      qtl_cs_variants_hg38
    ),
    
    chrombpnet_variants_supporting_qtl_cs_hg38 =
      collapse_variant_strings(
        chrombpnet_variants_qtl_cs_hg38[
          QTL_CS_has_ChromBPNet_variant
        ]
      ),
    
    n_chrombpnet_variants_supporting_qtl_cs =
      count_variants_in_field(
        collapse_variant_strings(
          chrombpnet_variants_qtl_cs_hg38[
            QTL_CS_has_ChromBPNet_variant
          ]
        )
      ),
    
    ## ---------------------------------------------------------------
    ## Colocalization details
    ## ---------------------------------------------------------------
    
    coloc_qtl_hit_variants_hg38 = collapse_variant_strings(
      coloc_qtl_hit_variants_hg38
    ),
    
    max_coloc_PP.H4.abf = max_or_na(
      max_coloc_PP.H4.abf
    ),
    
    coloc_supporting_gwas_studies = collapse_semicolon_fields(
      coloc_supporting_gwas_studies
    ),
    
    n_coloc_supporting_studies = count_semicolon_values(
      coloc_supporting_gwas_studies
    ),
    
    ## ---------------------------------------------------------------
    ## Match status
    ## ---------------------------------------------------------------
    
    all_qtl_cs_details_matched = all(
      qtl_cs_detail_matched
    ),
    
    n_qtl_cs_without_detail = sum(
      !qtl_cs_detail_matched
    ),
    
    .groups = "drop"
  )


## QC: the collapsed object must contain one row per unique trait-locus.
duplicate_qtl_summary_loci <- qtl_pair_summary_at_locus %>%
  count(
    across(all_of(KEY_LOCUS)),
    name = "n"
  ) %>%
  filter(n > 1)

stopifnot(nrow(duplicate_qtl_summary_loci) == 0)


## ======================================================================
## 14) SUPPLEMENTARY TABLE 8: ONE ROW PER UNIQUE SUPPORTED GWAS LOCUS
## ======================================================================
##
## Start from supported_coloc_loci, which is the exact unique-locus
## population used for the ChromBPNet + QTL colocalization portion of Plot 3.
##
## Do not expand supported_coloc_loci into individual QTL-pair rows.
## Instead, join the already-collapsed locus-level QTL summary.

supplementary_table8_cbpnet_molqtl <- supported_coloc_loci %>%
  left_join(
    gwas_cs_detail_trait_locus,
    by = KEY_LOCUS
  ) %>%
  left_join(
    qtl_pair_summary_at_locus,
    by = KEY_LOCUS
  ) %>%
  mutate(
    ## Combined ChromBPNet-supported variants from the GWAS CS and
    ## genuinely supporting paired QTL credible sets.
    chrombpnet_variants_all_supporting_cs_at_locus_hg38 =
      combine_variant_fields(
        chrombpnet_variants_gwas_cs_hg38,
        chrombpnet_variants_supporting_qtl_cs_hg38
      ),
    
    n_chrombpnet_variants_all_supporting_cs_at_locus =
      count_variants_in_field(
        chrombpnet_variants_all_supporting_cs_at_locus_hg38
      ),
    
    ## Readable locus-level source of ChromBPNet support.
    chrombpnet_support_source = case_when(
      locus_has_cbpnet_gwas_cs &
        n_chrombpnet_supporting_qtl_credible_sets > 0L ~
        "GWAS CS and paired QTL CS",
      
      locus_has_cbpnet_gwas_cs &
        n_chrombpnet_supporting_qtl_credible_sets == 0L ~
        "GWAS CS only",
      
      !locus_has_cbpnet_gwas_cs &
        n_chrombpnet_supporting_qtl_credible_sets > 0L ~
        "Paired QTL CS only",
      
      TRUE ~
        "No ChromBPNet support"
    )
  ) %>%
  rename(
    ## All GWAS datasets contributing the collapsed trait-locus.
    GWAS_studies_at_locus =
      contributing_gwas_studies,
    
    ## GWAS datasets in which a retained significant coloc pairing occurred.
    GWAS_studies_supporting_colocalization =
      coloc_supporting_gwas_studies,
    
    n_GWAS_studies_at_locus =
      n_gwas_studies_at_locus,
    
    n_GWAS_studies_supporting_colocalization =
      n_coloc_supporting_studies
  ) %>%
  select(
    ## ---------------------------------------------------------------
    ## GWAS locus identifiers
    ## ---------------------------------------------------------------
    
    trait,
    region,
    cs_index,
    region_cs,
    
    ## ---------------------------------------------------------------
    ## GWAS study information
    ## ---------------------------------------------------------------
    
    GWAS_studies_at_locus,
    GWAS_studies_supporting_colocalization,
    n_GWAS_studies_at_locus,
    n_GWAS_studies_supporting_colocalization,
    
    ## ---------------------------------------------------------------
    ## Locus classification
    ## ---------------------------------------------------------------
    
    qtl_coloc_category,
    chrombpnet_support_source,
    
    ## ---------------------------------------------------------------
    ## Linked molecular features
    ## ---------------------------------------------------------------
    
    linked_eGenes,
    linked_caPeaks,
    n_linked_eGenes,
    n_linked_caPeaks,
    
    ## ---------------------------------------------------------------
    ## All paired QTL credible sets
    ## ---------------------------------------------------------------
    
    all_paired_qtl_credible_sets,
    paired_eQTL_credible_sets,
    paired_caQTL_credible_sets,
    n_paired_qtl_credible_sets,
    
    ## ---------------------------------------------------------------
    ## QTL CSs that actually provide ChromBPNet support
    ## ---------------------------------------------------------------
    
    chrombpnet_supporting_qtl_credible_sets,
    chrombpnet_supporting_eQTL_credible_sets,
    chrombpnet_supporting_caQTL_credible_sets,
    n_chrombpnet_supporting_qtl_credible_sets,
    
    ## ---------------------------------------------------------------
    ## GWAS fine-mapping details
    ## ---------------------------------------------------------------
    
    gwas_lead_variants_hg37,
    gwas_lead_variants_by_study_hg37,
    maximum_gwas_lead_pip,
    gwas_cs_variants_hg37,
    
    ## ---------------------------------------------------------------
    ## QTL fine-mapping and coloc details
    ## ---------------------------------------------------------------
    
    qtl_lead_variants_hg38,
    maximum_qtl_lead_pip,
    coloc_qtl_hit_variants_hg38,
    max_coloc_PP.H4.abf,
    all_paired_qtl_cs_variants_hg38,
    
    ## ---------------------------------------------------------------
    ## ChromBPNet-supported variants
    ## ---------------------------------------------------------------
    
    chrombpnet_variants_gwas_cs_hg37,
    chrombpnet_variants_gwas_cs_hg38,
    chrombpnet_variants_supporting_qtl_cs_hg38,
    chrombpnet_variants_all_supporting_cs_at_locus_hg38,
    
    n_chrombpnet_variants_gwas_cs,
    n_chrombpnet_variants_supporting_qtl_cs,
    n_chrombpnet_variants_all_supporting_cs_at_locus,
    
    ## ---------------------------------------------------------------
    ## QC columns
    ## ---------------------------------------------------------------
    
    all_qtl_cs_details_matched,
    n_qtl_cs_without_detail
  ) %>%
  arrange(
    trait,
    region,
    cs_index
  )


## ----------------------------------------------------------------------
## Compact publication-facing version
## ----------------------------------------------------------------------

supplementary_table8_compact <- supplementary_table8_cbpnet_molqtl %>%
  select(
    trait,
    region_cs,
    
    ## GWAS study information
    GWAS_studies_at_locus,
    GWAS_studies_supporting_colocalization,
    
    ## Locus classification
    qtl_coloc_category,
    chrombpnet_support_source,
    
    ## Molecular features
    linked_eGenes,
    linked_caPeaks,
    
    ## QTL credible sets
    all_paired_qtl_credible_sets,
    chrombpnet_supporting_qtl_credible_sets,
    
    ## Fine-mapping and coloc details
    gwas_lead_variants_hg37,
    qtl_lead_variants_hg38,
    coloc_qtl_hit_variants_hg38,
    max_coloc_PP.H4.abf,
    
    ## ChromBPNet-supported variants
    chrombpnet_variants_gwas_cs_hg37,
    chrombpnet_variants_gwas_cs_hg38,
    chrombpnet_variants_supporting_qtl_cs_hg38,
    chrombpnet_variants_all_supporting_cs_at_locus_hg38
  )


## ======================================================================
## 15) FINAL QC: NO DUPLICATED LOCI AND EXACT AGREEMENT WITH PLOT 3
## ======================================================================

## ----------------------------------------------------------------------
## 15A) Explicit duplicate check
## ----------------------------------------------------------------------

supp_table_duplicate_qc <- supplementary_table8_cbpnet_molqtl %>%
  count(
    trait,
    region,
    cs_index,
    region_cs,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1)

write.csv(
  supp_table_duplicate_qc,
  file.path(
    output_dir,
    "QC_duplicate_supplementary_table8_loci.csv"
  ),
  row.names = FALSE
)

if (nrow(supp_table_duplicate_qc) > 0) {
  stop(
    "Supplementary Table 8 contains duplicated trait-GWAS loci. ",
    "Inspect QC_duplicate_supplementary_table8_loci.csv."
  )
}


## ----------------------------------------------------------------------
## 15B) Table rows must exactly equal supported plot loci
## ----------------------------------------------------------------------

table_locus_keys <- supplementary_table8_cbpnet_molqtl %>%
  select(all_of(KEY_LOCUS)) %>%
  distinct()

plot_locus_keys <- supported_coloc_loci %>%
  select(all_of(KEY_LOCUS)) %>%
  distinct()


extra_supported_loci_in_table <- table_locus_keys %>%
  anti_join(
    plot_locus_keys,
    by = KEY_LOCUS
  )


missing_supported_loci_from_table <- plot_locus_keys %>%
  anti_join(
    table_locus_keys,
    by = KEY_LOCUS
  )


write.csv(
  extra_supported_loci_in_table,
  file.path(
    output_dir,
    "QC_extra_supported_loci_in_supp_table8.csv"
  ),
  row.names = FALSE
)

write.csv(
  missing_supported_loci_from_table,
  file.path(
    output_dir,
    "QC_missing_supported_loci_from_supp_table8.csv"
  ),
  row.names = FALSE
)


if (
  nrow(extra_supported_loci_in_table) > 0 ||
  nrow(missing_supported_loci_from_table) > 0
) {
  stop(
    "Supplementary Table 8 locus membership does not exactly match ",
    "supported_coloc_loci. Inspect the extra/missing locus QC files."
  )
}


## The row count must equal the number of unique supported plot loci.
stopifnot(
  nrow(supplementary_table8_cbpnet_molqtl) ==
    nrow(plot_locus_keys)
)


## Every table row must represent a unique locus.
stopifnot(
  nrow(supplementary_table8_cbpnet_molqtl) ==
    nrow(table_locus_keys)
)


## ----------------------------------------------------------------------
## 15C) Trait-specific table counts must reproduce Plot 3
## ----------------------------------------------------------------------

supp_locus_qc <- supplementary_table8_cbpnet_molqtl %>%
  count(
    trait,
    name = "n_supported_coloc_loci_in_supp_table"
  ) %>%
  full_join(
    trait_summary %>%
      select(
        trait,
        expected_supported_coloc_loci =
          cs_qtl_coloc_with_cbpnet
      ),
    by = "trait"
  ) %>%
  mutate(
    n_supported_coloc_loci_in_supp_table = replace_na(
      n_supported_coloc_loci_in_supp_table,
      0L
    ),
    
    expected_supported_coloc_loci = replace_na(
      expected_supported_coloc_loci,
      0L
    ),
    
    difference =
      n_supported_coloc_loci_in_supp_table -
      expected_supported_coloc_loci,
    
    passes_check = difference == 0
  )


write.csv(
  supp_locus_qc,
  file.path(
    output_dir,
    "QC_supp_table8_vs_plot3_locus_counts.csv"
  ),
  row.names = FALSE
)


if (any(!supp_locus_qc$passes_check)) {
  stop(
    "Supplementary Table 8 does not reproduce the supported-colocalized ",
    "locus counts used in Plot 3. Inspect ",
    "QC_supp_table8_vs_plot3_locus_counts.csv."
  )
}


## ----------------------------------------------------------------------
## 15D) Confirm support-source labels agree with actual variant fields
## ----------------------------------------------------------------------

supp_support_source_qc <- supplementary_table8_cbpnet_molqtl %>%
  mutate(
    actual_gwas_support =
      !is.na(chrombpnet_variants_gwas_cs_hg38) &
      chrombpnet_variants_gwas_cs_hg38 != "",
    
    actual_qtl_support =
      !is.na(chrombpnet_variants_supporting_qtl_cs_hg38) &
      chrombpnet_variants_supporting_qtl_cs_hg38 != "",
    
    expected_support_source = case_when(
      actual_gwas_support & actual_qtl_support ~
        "GWAS CS and paired QTL CS",
      
      actual_gwas_support & !actual_qtl_support ~
        "GWAS CS only",
      
      !actual_gwas_support & actual_qtl_support ~
        "Paired QTL CS only",
      
      TRUE ~
        "No ChromBPNet support"
    ),
    
    passes_check =
      chrombpnet_support_source ==
      expected_support_source
  )


write.csv(
  supp_support_source_qc %>%
    filter(!passes_check),
  file.path(
    output_dir,
    "QC_supp_table8_support_source_disagreement.csv"
  ),
  row.names = FALSE
)


if (any(!supp_support_source_qc$passes_check)) {
  stop(
    "One or more Supplementary Table 8 support-source labels disagree ",
    "with the actual ChromBPNet variant fields. Inspect ",
    "QC_supp_table8_support_source_disagreement.csv."
  )
}


## ----------------------------------------------------------------------
## 15E) Confirm supporting GWAS study information is present
## ----------------------------------------------------------------------

missing_supporting_gwas_study_qc <-
  supplementary_table8_cbpnet_molqtl %>%
  filter(
    is.na(GWAS_studies_supporting_colocalization) |
      GWAS_studies_supporting_colocalization == ""
  )


write.csv(
  missing_supporting_gwas_study_qc,
  file.path(
    output_dir,
    "QC_missing_supporting_GWAS_studies.csv"
  ),
  row.names = FALSE
)


if (nrow(missing_supporting_gwas_study_qc) > 0) {
  stop(
    nrow(missing_supporting_gwas_study_qc),
    " retained loci do not have a GWAS study supporting colocalization. ",
    "Inspect QC_missing_supporting_GWAS_studies.csv."
  )
}


## ----------------------------------------------------------------------
## 15F) Supporting studies must be included among studies at the locus
## ----------------------------------------------------------------------

supporting_study_membership_qc <-
  supplementary_table8_cbpnet_molqtl %>%
  rowwise() %>%
  mutate(
    studies_at_locus_list = list(
      if (
        is.na(GWAS_studies_at_locus) ||
        GWAS_studies_at_locus == ""
      ) {
        character(0)
      } else {
        trimws(
          unlist(
            strsplit(
              GWAS_studies_at_locus,
              ";",
              fixed = TRUE
            )
          )
        )
      }
    ),
    
    supporting_studies_list = list(
      if (
        is.na(GWAS_studies_supporting_colocalization) ||
        GWAS_studies_supporting_colocalization == ""
      ) {
        character(0)
      } else {
        trimws(
          unlist(
            strsplit(
              GWAS_studies_supporting_colocalization,
              ";",
              fixed = TRUE
            )
          )
        )
      }
    ),
    
    supporting_studies_are_at_locus =
      all(
        supporting_studies_list %in%
          studies_at_locus_list
      )
  ) %>%
  ungroup()


write.csv(
  supporting_study_membership_qc %>%
    filter(!supporting_studies_are_at_locus) %>%
    select(
      trait,
      region,
      cs_index,
      region_cs,
      GWAS_studies_at_locus,
      GWAS_studies_supporting_colocalization
    ),
  file.path(
    output_dir,
    "QC_supporting_GWAS_study_not_at_locus.csv"
  ),
  row.names = FALSE
)


if (
  any(
    !supporting_study_membership_qc$
    supporting_studies_are_at_locus
  )
) {
  stop(
    "At least one GWAS study supporting colocalization is not included ",
    "among the GWAS studies assigned to the locus. Inspect ",
    "QC_supporting_GWAS_study_not_at_locus.csv."
  )
}


## ----------------------------------------------------------------------
## 15G) Report missing QTL-CS detail without changing locus counts
## ----------------------------------------------------------------------

missing_qtl_detail_in_supp_table <-
  supplementary_table8_cbpnet_molqtl %>%
  filter(
    !all_qtl_cs_details_matched |
      n_qtl_cs_without_detail > 0
  )


write.csv(
  missing_qtl_detail_in_supp_table,
  file.path(
    output_dir,
    "QC_missing_qtl_detail_in_supp_table.csv"
  ),
  row.names = FALSE
)


if (nrow(missing_qtl_detail_in_supp_table) > 0) {
  warning(
    nrow(missing_qtl_detail_in_supp_table),
    " supported loci contain at least one paired QTL CS without matched ",
    "descriptive details. Locus-level counts remain valid, but inspect ",
    "QC_missing_qtl_detail_in_supp_table.csv before publication."
  )
}


## ======================================================================
## 16) EXPORT ANALYSIS TABLES AND SUPPLEMENTARY TABLE 8
## ======================================================================

new_percent_by_trait <- trait_summary %>%
  transmute(
    trait,
    trait_label,
    
    total_unique_gwas_loci =
      total_cs_all,
    
    loci_with_cbpnet_support =
      cs_with_cbpnet_support,
    
    loci_with_qtl_coloc_and_cbpnet_support =
      cs_qtl_coloc_with_cbpnet,
    
    percent_gwas_loci_with_cbpnet_support =
      ifelse(
        total_cs_all > 0,
        100 *
          cs_with_cbpnet_support /
          total_cs_all,
        NA_real_
      ),
    
    percent_gwas_loci_with_qtl_coloc_and_cbpnet_support =
      ifelse(
        total_cs_all > 0,
        100 *
          cs_qtl_coloc_with_cbpnet /
          total_cs_all,
        NA_real_
      )
  )


new_percent_summary_all_traits <- new_percent_by_trait %>%
  summarise(
    trait_set = "All traits",
    
    n_traits = n(),
    
    mean_percent_gwas_loci_with_cbpnet_support =
      mean(
        percent_gwas_loci_with_cbpnet_support,
        na.rm = TRUE
      ),
    
    mean_percent_gwas_loci_with_qtl_coloc_and_cbpnet_support =
      mean(
        percent_gwas_loci_with_qtl_coloc_and_cbpnet_support,
        na.rm = TRUE
      ),
    
    pooled_percent_gwas_loci_with_cbpnet_support =
      ifelse(
        sum(total_unique_gwas_loci) > 0,
        100 *
          sum(loci_with_cbpnet_support) /
          sum(total_unique_gwas_loci),
        NA_real_
      ),
    
    pooled_percent_gwas_loci_with_qtl_coloc_and_cbpnet_support =
      ifelse(
        sum(total_unique_gwas_loci) > 0,
        100 *
          sum(
            loci_with_qtl_coloc_and_cbpnet_support
          ) /
          sum(total_unique_gwas_loci),
        NA_real_
      )
  )


new_percent_summary_coloc_traits <- new_percent_by_trait %>%
  filter(trait %in% keep_traits) %>%
  summarise(
    trait_set =
      "Traits with at least one QTL-GWAS colocalized locus",
    
    n_traits = n(),
    
    mean_percent_gwas_loci_with_cbpnet_support =
      mean(
        percent_gwas_loci_with_cbpnet_support,
        na.rm = TRUE
      ),
    
    mean_percent_gwas_loci_with_qtl_coloc_and_cbpnet_support =
      mean(
        percent_gwas_loci_with_qtl_coloc_and_cbpnet_support,
        na.rm = TRUE
      ),
    
    pooled_percent_gwas_loci_with_cbpnet_support =
      ifelse(
        sum(total_unique_gwas_loci) > 0,
        100 *
          sum(loci_with_cbpnet_support) /
          sum(total_unique_gwas_loci),
        NA_real_
      ),
    
    pooled_percent_gwas_loci_with_qtl_coloc_and_cbpnet_support =
      ifelse(
        sum(total_unique_gwas_loci) > 0,
        100 *
          sum(
            loci_with_qtl_coloc_and_cbpnet_support
          ) /
          sum(total_unique_gwas_loci),
        NA_real_
      )
  )


new_percent_summary <- bind_rows(
  new_percent_summary_all_traits,
  new_percent_summary_coloc_traits
)




## ======================================================================
## SUPPLEMENTARY TABLE: CHROMBPNET-SUPPORTED GWAS LOCI WITHOUT QTL COLOC
## ======================================================================
##
## This table corresponds exactly to the "ChromBPNet only" category in
## Plot 3:
##
##   cbpnet_category == "ChromBPNet_specific"
##
## Each row represents one unique trait-GWAS credible-set locus.
## These loci contain at least one ChromBPNet-prioritized variant in the
## GWAS credible set but have no retained eQTL-GWAS or caQTL-GWAS
## colocalization.

chrombpnet_only_gwas_loci <- trait_locus_master %>%
  filter(
    cbpnet_category == "ChromBPNet_specific"
  ) %>%
  select(
    trait,
    region,
    cs_index,
    region_cs,
    n_studies_at_locus,
    study_ids,
    has_cbpnet_gwas_cs,
    qtl_coloc_cs,
    cbpnet_category
  ) %>%
  left_join(
    gwas_cs_detail_trait_locus,
    by = KEY_LOCUS
  ) %>%
  mutate(
    chrombpnet_support_source = "GWAS CS only"
  ) %>%
  rename(
    GWAS_studies_at_locus = contributing_gwas_studies,
    n_GWAS_studies_at_locus = n_gwas_studies_at_locus
  ) %>%
  select(
    ## GWAS locus identifiers
    trait,
    region,
    cs_index,
    region_cs,
    
    ## GWAS study information
    GWAS_studies_at_locus,
    n_GWAS_studies_at_locus,
    
    ## ChromBPNet classification
    chrombpnet_support_source,
    
    ## GWAS fine-mapping details
    gwas_lead_variants_hg37,
    gwas_lead_variants_by_study_hg37,
    maximum_gwas_lead_pip,
    gwas_cs_variants_hg37,
    
    ## ChromBPNet-supported variants
    chrombpnet_variants_gwas_cs_hg37,
    chrombpnet_variants_gwas_cs_hg38,
    n_chrombpnet_variants_gwas_cs
  ) %>%
  arrange(
    trait,
    region,
    cs_index
  )

## ----------------------------------------------------------------------
## QC: exactly one row per unique trait-GWAS locus
## ----------------------------------------------------------------------

chrombpnet_only_duplicate_qc <- chrombpnet_only_gwas_loci %>%
  count(
    trait,
    region,
    cs_index,
    region_cs,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1)

if (nrow(chrombpnet_only_duplicate_qc) > 0) {
  write.csv(
    chrombpnet_only_duplicate_qc,
    file.path(
      output_dir,
      "QC_duplicate_ChromBPNet_only_GWAS_loci.csv"
    ),
    row.names = FALSE
  )
  
  stop(
    "The ChromBPNet-only GWAS table contains duplicated trait-locus rows. ",
    "Inspect QC_duplicate_ChromBPNet_only_GWAS_loci.csv."
  )
}

chrombpnet_only_plot_qc <- chrombpnet_only_gwas_loci %>%
  count(
    trait,
    name = "n_chrombpnet_only_loci_in_table"
  ) %>%
  full_join(
    trait_summary %>%
      select(
        trait,
        expected_chrombpnet_only_loci = cs_cbpnet_specific
      ),
    by = "trait"
  ) %>%
  mutate(
    n_chrombpnet_only_loci_in_table = replace_na(
      n_chrombpnet_only_loci_in_table,
      0L
    ),
    expected_chrombpnet_only_loci = replace_na(
      expected_chrombpnet_only_loci,
      0L
    ),
    difference =
      n_chrombpnet_only_loci_in_table -
      expected_chrombpnet_only_loci,
    passes_check = difference == 0
  )

write.csv(
  chrombpnet_only_plot_qc,
  file.path(
    output_dir,
    "QC_ChromBPNet_only_GWAS_table_vs_plot3.csv"
  ),
  row.names = FALSE
)

if (any(!chrombpnet_only_plot_qc$passes_check)) {
  stop(
    "The ChromBPNet-only GWAS table does not reproduce the blue-category ",
    "counts in Plot 3. Inspect QC_ChromBPNet_only_GWAS_table_vs_plot3.csv."
  )
}

chrombpnet_only_gwas_loci_compact <- chrombpnet_only_gwas_loci %>%
  select(
    trait,
    GWAS_studies_at_locus,
    region_cs,
    gwas_lead_variants_hg37,
    maximum_gwas_lead_pip,
    chrombpnet_variants_gwas_cs_hg37,
    chrombpnet_variants_gwas_cs_hg38,
    n_chrombpnet_variants_gwas_cs
  )

chrombpnet_only_compact_file <- file.path(
  output_dir,
  "Supplementary_Table_ChromBPNet_only_GWAS_loci_compact.csv"
)

write.csv(
  chrombpnet_only_gwas_loci_compact,
  chrombpnet_only_compact_file,
  row.names = FALSE
)

## ----------------------------------------------------------------------
## Export primary analysis objects
## ----------------------------------------------------------------------

write.csv(
  trait_locus_master,
  file.path(
    output_dir,
    "trait_locus_master_for_plots_and_supp_table8.csv"
  ),
  row.names = FALSE
)

write.csv(
  coloc_table_summary,
  file.path(
    output_dir,
    "plot2_qtl_gwas_coloc_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  trait_summary,
  file.path(
    output_dir,
    "plot3_chrombpnet_support_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  new_percent_by_trait,
  file.path(
    output_dir,
    "chrombpnet_percent_by_trait.csv"
  ),
  row.names = FALSE
)

write.csv(
  new_percent_summary,
  file.path(
    output_dir,
    "chrombpnet_percent_summary.csv"
  ),
  row.names = FALSE
)


## ----------------------------------------------------------------------
## Export locus-level QTL pairing summary
## ----------------------------------------------------------------------

write.csv(
  qtl_pair_summary_at_locus,
  file.path(
    output_dir,
    "Supplementary_Table_8_QTL_pair_summary_by_locus.csv"
  ),
  row.names = FALSE
)


## ----------------------------------------------------------------------
## Export Supplementary Table 8
## ----------------------------------------------------------------------

full_supp_table_file <- file.path(
  output_dir,
  paste0(
    "Supplementary_Table_8_",
    "ChromBPNet_supported_QTL_GWAS_loci_full.csv"
  )
)

compact_supp_table_file <- file.path(
  output_dir,
  paste0(
    "Supplementary_Table_8_",
    "ChromBPNet_supported_QTL_GWAS_loci_compact.csv"
  )
)


write.csv(
  supplementary_table8_cbpnet_molqtl,
  full_supp_table_file,
  row.names = FALSE
)

write.csv(
  supplementary_table8_compact,
  compact_supp_table_file,
  row.names = FALSE
)


message("Unified analysis completed successfully.")

message(
  "Unique supported loci in Plot 3: ",
  nrow(supported_coloc_loci)
)

message(
  "Rows in full Supplementary Table 8: ",
  nrow(supplementary_table8_cbpnet_molqtl)
)

message(
  "Rows in compact Supplementary Table 8: ",
  nrow(supplementary_table8_compact)
)

message(
  "Loci missing a supporting GWAS study: ",
  nrow(missing_supporting_gwas_study_qc)
)

message(
  "Master locus table: ",
  file.path(
    output_dir,
    "trait_locus_master_for_plots_and_supp_table8.csv"
  )
)

message(
  "Full Supplementary Table 8: ",
  full_supp_table_file
)

message(
  "Compact Supplementary Table 8: ",
  compact_supp_table_file
)

message(
  "Combined plot: ",
  plot_output_file
)