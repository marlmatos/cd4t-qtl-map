rm(list = ls())
options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(ggplot2)

setwd("/gchm_lab/collab/marlis_pj/coloc")

files = system("ls coloc_results/eqtl_gwas_coloc/*.csv", intern = T)
basenames = basename(files)

# Extract the year (start of filename)
years <- str_extract(basenames, "[0-9]{4}")

# Extract the TRAIT_POP (last two parts before "_preprocessed")
traits <- str_extract(basenames, "[^_/]+_[^_/]+(?=_preprocessed)")

# Combine
combined <- paste(years, traits, sep = "_")
print(sort(combined))

summarise_coloc = function(x){
  
  print(files[x])
  
  gwas_coloc = fread(files[x]) %>%
    filter(pval < 5e-08 & pval_nominal < 1e-03)
  
  out_df = data.frame(gwas = combined[x],
                      n_eqtls = nrow(gwas_coloc), # No coloc eQTLs
                      n_genes = length(unique(gwas_coloc$gene))) # No coloc genes
  
  return(out_df)
}

results = lapply(1:length(files), function(x) summarise_coloc(x))

summary_df = dplyr::bind_rows(results) %>% 
  arrange(desc(n_eqtls)) %>%
  distinct(gwas, .keep_all=T) # remove same gwas studies from same year

top_20 = summary_df[1:20,]
top_20$gwas <- factor(top_20$gwas, levels = rev(unique(top_20$gwas)))

# Make stacked barplot
# Barplot
top_long <- pivot_longer(
  top_20,
  cols = c(n_eqtls, n_genes),
  names_to = "type",
  values_to = "count"
)

png("plots/gwas_eqtl_top20_barplot.png", width = 8, height = 6, units = "in", res=300)

ggplot(top_long, aes(x = gwas, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "GWAS", y = "Count", fill = "Measure") +
  theme_minimal() + 
  labs(x = "GWAS", y = "Count", fill = "") # ← rename legend title here

dev.off()

# Same for caQTL
files = system("ls coloc_results/ca_gwas_coloc/*.csv", intern = T)
basenames = basename(files)
basenames = sub("CD4T_chromatin_", "", basenames)

# Extract the year (start of filename)
years <- str_extract(basenames, "[0-9]{4}")

# Extract the TRAIT_POP (last two parts before "_preprocessed")
traits <- str_extract(basenames, "[^_/]+_[^_/]+(?=_preprocessed)")

# Combine
combined <- paste(years, traits, sep = "_")
print(sort(combined))

summarise_coloc = function(x){
  
  print(files[x])
  
  gwas_coloc = fread(files[x]) %>%
    filter(pval < 5e-08 & pval_nominal < 1e-03)
  
  out_df = data.frame(gwas = combined[x],
                      n_caqtls = nrow(gwas_coloc), # No coloc eQTLs
                      n_peaks = length(unique(gwas_coloc$peak))) # No coloc genes
  
  return(out_df)
}

results2 = lapply(1:length(files), function(x) summarise_coloc(x))

summary_df = dplyr::bind_rows(results2) %>% 
  arrange(desc(n_caqtls)) %>%
  distinct(gwas, .keep_all=T) # remove same gwas studies from same year

top_20 = summary_df[1:20,]
top_20$gwas <- factor(top_20$gwas, levels = rev(unique(top_20$gwas)))

# Make stacked barplot
# Barplot
top_long <- pivot_longer(
  top_20,
  cols = c(n_caqtls, n_peaks),
  names_to = "type",
  values_to = "count"
)

png("plots/gwas_caqtl_top20_barplot.png", width = 8, height = 6, units = "in", res=300)

ggplot(top_long, aes(x = gwas, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "GWAS", y = "Count", fill = "Measure") +
  theme_minimal() + 
  labs(x = "GWAS", y = "Count", fill = "") # ← rename legend title here

dev.off()

