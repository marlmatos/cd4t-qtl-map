rm(list=ls())
options(bitmapType = "cairo")
library(ggplot2)
library(data.table)
library(tidyverse)

setwd("/gchm_lab/collab/marlis_pj/coloc")

# folder names
file_names = unlist(fread("preprocessed_folder_names.txt", header = F))

gwas_id = file_names[58]

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
  
  gwas_dir <- '/gchm_lab/collab_gwas/finemapping_autoimmune/'
  
  # Get region RDS files
  regions <- list.files(
    path = file.path(gwas_dir, "output/nathan_completed", gwas_id),
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  regions <- regions[grepl("_1.5Mb/susie/", regions)]
  gwas_cs <- sapply(regions, count_cs)
  total_cs <- sum(gwas_cs)
  
  # Check if coloc result files exist
  eqtl_path  <- paste0("coloc_results/eqtl_gwas_coloc/", gwas_id, "_coloc_results.csv")
  caqtl_path <- paste0("coloc_results/ca_gwas_coloc/CD4T_chromatin_", gwas_id, "_coloc_results.csv")
  
  if (!file.exists(eqtl_path) || !file.exists(caqtl_path)) return(NULL)
  
  # Read coloc results
  eqtl_coloc_result <- fread(eqtl_path) %>%
    mutate(region_cs = paste(region, idx1, sep = "_")) %>%
    filter(pval < 1e-5, pval_nominal < 1e-3)
  
  caqtl_coloc_result <- fread(caqtl_path) %>%
    mutate(region_cs = paste(region, idx1, sep = "_")) %>%
    filter(pval < 1e-5, pval_nominal < 1e-3)
  
  shared_coloc_result <- caqtl_coloc_result %>%
    filter(region_cs %in% eqtl_coloc_result$region_cs) %>%
    mutate(gwas_id = gwas_id)
  
  # Summary stats
  only_eqtl <- length(setdiff(unique(eqtl_coloc_result$region_cs), caqtl_coloc_result$region_cs))
  only_caqtl <- length(setdiff(unique(caqtl_coloc_result$region_cs), eqtl_coloc_result$region_cs))
  both <- length(unique(shared_coloc_result$region_cs))
  
  summary_table <- data.frame(
    trait = gwas_id,
    total_cond_indep = total_cs,
    just_eqtl = only_eqtl,
    just_caqtl = only_caqtl,
    both = both
  )
  
  return(list(summary = summary_table, shared = shared_coloc_result))
}


combined_table = lapply(1:length(file_names), function(x) make_table(file_names[x]))

results_filtered <- combined_table[!sapply(combined_table, function(x) {
  is.null(x) || (is.atomic(x) && any(is.na(x))) || (is.data.frame(x) && ncol(x) == 0)
})]

# Combine and save summary table
results2 <- bind_rows(lapply(results_filtered, `[[`, "summary"))
fwrite(results2, "~/cd4_qtl_paper_figures/figure_4/data/CD4T_coloc_summary_table.csv")

# Combine and save shared coloc results
all_shared <- bind_rows(lapply(results_filtered, `[[`, "shared"))
fwrite(all_shared, "~/cd4_qtl_paper_figures/figure_4/data/CD4T_all_shared_coloc.csv")


# # make plot
# coloc_table = fread("coloc_results_barplot_df.txt") %>%
#   arrange(desc(both))
# 
# # take first 10
# coloc_table = coloc_table[1:20,]
# 
# plot_df = pivot_longer(
#   coloc_table,
#   cols = c(just_eqtl, just_caqtl, both),
#   names_to = "category",
#   values_to = "count"
# )
# plot_df$proportion = plot_df$count / plot_df$total_cond_indep
# 
# library(ggplot2)
# 
# ggplot(plot_df, aes(x=trait, y=proportion, fill=category)) +
#   geom_bar(stat="identity") +
#   geom_text(
#     aes(label=count),
#     position=position_stack(vjust=0.5), # center text inside the bar segment
#     color="white",
#     size=4
#   ) +
#   labs(
#     title="Proportion of Conditional Independent Variants",
#     y="Proportion",
#     x="Trait"
#   ) +
#   theme_cowplot() +
#   scale_y_continuous(labels=scales::percent_format()) +
#   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# 
