.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))

options(bitmapType="cairo")
library(data.table)
library(tidyverse)
library(ggplot2)
library(readr)
library(dplyr)

metadata<-read_delim("~/cd4_qtl_paper_figures/figure_6/gwas_traits_stingseq.csv")
setwd("~/cd4_qtl_paper_figures/figure_6")

files_caqtl = system("ls /gcgl/sghatan/marlis_pj/coloc/coloc_results/ca_gwas_coloc/*.csv", intern = T)
files_eqtl = system("ls /gcgl/sghatan/marlis_pj/coloc/coloc_results/eqtl_gwas_coloc/*.csv", intern = T)


# From: CD4T_chromatin_<YEAR>_<ID>_<TRAIT>_<POP>_preprocessed...
year_caqtl   <- str_extract(files_caqtl, "(?<=CD4T_chromatin_)[0-9]{4}")
trait_caqtl  <- str_extract(files_caqtl, "[^_/]+_[^_]+_[^_]+(?=_preprocessed)")  # RA_EUR, etc.
gwas_caqtl   <- paste(year_caqtl, trait_caqtl, sep = "_")

# From: 2024_PanUKBB_RA_EUR_preprocessed_coloc_results.csv
year_eqtl   <- str_extract(files_eqtl, "(?<=/)[0-9]{4}")
trait_eqtl  <- str_extract(files_eqtl, "[^_/]+_[^_/]+_[^_/]+(?=_preprocessed)")  # PanUKBB_RA_EUR, etc.
gwas_eqtl   <- paste(year_eqtl, trait_eqtl, sep = "_")

##caQTL_gwas_coloc_res
summarise_caqtl <- function(x) {
  gwas_coloc <- fread(files_caqtl[x]) %>%
    filter(pval < 5e-08 & pval_nominal < 1e-03) %>%
    arrange(desc(PP.H4.abf)) %>%
    mutate(unique_id = paste(GWAS_variant_GRC37, eQTL_variant_GRC37, peak, sep = "_")) %>%
    distinct(unique_id, .keep_all = TRUE)
  
  data.frame(
    gwas = gwas_caqtl[x],
    n_caqtls = nrow(gwas_coloc),
    n_peaks = length(unique(gwas_coloc$peak))
  )
}

##eQTL_gwas_coloc_res
summarise_eqtl <- function(x) {
  gwas_coloc <- fread(files_eqtl[x]) %>%
    filter(pval < 5e-08 & pval_nominal < 1e-03) %>%
    arrange(desc(PP.H4.abf)) %>%
    mutate(unique_id = paste(GWAS_variant_GRC37, eQTL_variant_GRC37, gene, sep = "_")) %>%
    distinct(unique_id, .keep_all = TRUE)
  
  data.frame(
    gwas = gwas_eqtl[x],
    n_eqtls = nrow(gwas_coloc),
    n_genes = length(unique(gwas_coloc$gene))
  )
}

caqtl_summary <- bind_rows(lapply(seq_along(files_caqtl), summarise_caqtl))
eqtl_summary  <- bind_rows(lapply(seq_along(files_eqtl), summarise_eqtl))
coloc_summary <- full_join(caqtl_summary, eqtl_summary, by = "gwas") %>%
  arrange(desc(n_eqtls)) %>%
  distinct(gwas, .keep_all=T) %>%  # remove same gwas studies from same year
  left_join(metadata, join_by(gwas==base_name))

top_20 = coloc_summary[1:20,]
top_20$gwas <- factor(top_20$gwas, levels = rev(unique(top_20$gwas)))

# Make stacked barplot
# Barplot
top_long <- pivot_longer(
  top_20,
  cols = c(n_eqtls, n_caqtls),
  names_to = "type",
  values_to = "count"
) 
top_long$type <- top_long$type %>%
  recode(n_eqtls = "eQTLs", n_caqtls = "caQTLs") %>%
  factor(levels = c("eQTLs", "caQTLs"))

missing_traits_names<-c(32296059 = "Asthma",
                        32581359= "Autoimmune thyroid disease",
                        36333501= "Rheumatoid Arthritis",
                        26394269="Primary biliary cirrhosis",
                        29083406="Asthma",
                        38296975="Systemic Sclerosis",
                        26482879="Atopic Dermatitis",
                        27723758="IgA Deficiency",
                        33493351="Systemic Lupus Erythematosus",
                        37156999="Inflamatory Bowel Disease"
                        
)


missing_traits_names2<-c("2024_PanUKBB_AST_EUR"="Asthma",
                         "2020_32296059_AST_EUR"="Asthma",
                         "2024_38982041_AITD_EUR"="Autoimmune thyroid disease")

top_long <- top_long %>% 
  left_join(lookup_tbl, by = "gwas") %>% 
  mutate(
    disease_phenotype = coalesce(disease_phenotype, disease_new_name)
  ) %>% 
  select(-disease_new_name)


# Color map
cmap <- c("caQTLs" = "#9CCB86",
          "eQTLs"  = "#C70E7B")

# Plot
png("~/cd4_qtl_paper_figures/figure_6gwas_eqtl_top20_barplot.pdf", width = 6, height = 5, units = "in", res=300)
ggplot(top_long, aes(x = gwas, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", color="gray20", size=0.2) +
  scale_fill_manual(values = cmap) +
  coord_flip() +
  theme_linedraw() + theme(panel.grid.minor.x = element_blank(),
                           panel.grid.minor.y = element_blank())+
  labs(x = "GWAS Trait", y = "Number of Colocalized Loci", fill = "")
dev.off()
