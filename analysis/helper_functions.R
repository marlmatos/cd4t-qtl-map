library(dplyr)
library(ggplot2)
library(readr) 
library(stringr)
#############################################################################
#prepare finemapping dataframes for integration
prep_finemap_df <- function(df, phenotype_col = NULL) {
  
  ## 1. choose the phenotype column -----------------------------------------
  if (is.null(phenotype_col)) {
    if ("gene" %in% names(df))      phenotype_col <- "gene"
    else if ("peak" %in% names(df)) phenotype_col <- "peak"
    else stop(
      "Can't determine the phenotype column. ",
      "Call with phenotype_col = 'gene' or 'peak'."
    )
  }
  if (!phenotype_col %in% names(df))
    stop("Column ‘", phenotype_col, "’ not found in the data frame.")
  
  ## 2. build / fix the required columns ------------------------------------
  df %>% 
    mutate(
      # numeric coordinate ----------------------------------------------------
      variant_pos = if ("variant_pos" %in% names(.))
        as.numeric(variant_pos)
      else as.numeric(sub(".*:(\\d+)\\[.*", "\\1", variant_id)),
      
      # chromosome string (ensure 'chr' prefix) -------------------------------
      chr = case_when(
        "chr" %in% names(.) ~ if_else(
          str_starts(chr, "chr"),
          as.character(chr),
          paste0("chr", chr)
        ),
        TRUE ~ paste0(
          "chr",
          sub("^([0-9XYMT]+):.*", "\\1", variant_id)
        )
      ),
      
      # credible-set identifier ----------------------------------------------
      finemapped_cs = if ("finemapped_cs" %in% names(.))
        finemapped_cs
      else paste(.data[[phenotype_col]], region, cs, sep = "_")
    )
}

################################################################################
####### Collapse each CS to an interval + variant list ────────────────────────
collapse_cs <- function(df) {
  df %>% 
    mutate(var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])", 
                         "chr\\1_\\2_\\3_\\4", 
                         variant_id)) %>%
    group_by(finemapped_cs, chr) %>% 
    summarise(
      start    = min(variant_pos),
      end      = max(variant_pos),
      variants = list(unique(variant_id)),
      .groups  = "drop"
    )
}