################################
## Script to Process overlap MosDisco Motif Instances
## Author: Marliette Matos
## Date: 10/13/2025
################################s
library(readr)
library(dplyr)
library(tidyr)
library(rlang)  
library(cowplot)
library(GenomicRanges)
library(igraph)
library(stringr)
library(tibble)
library(ggplot2)
options(bitmapType = "cairo")

variant_hits<-read_tsv("~/cd4_chrombpnet/chrombpnet_model_b7/motif_hit_calls/motifs_hits_FILTVARS/variant_hit_calls.tsv")
head(variant_hits)

#number of motif hits per allele
table(variant_hits$allele)

#how many distinct variants overlap a motif?
n_distinct(variant_hits$variant_loc)

#is the variant location always within the motif's start and end?
table((variant_hits$variant_loc >= variant_hits$start) & (variant_hits$variant_loc <= variant_hits$end)) 
#TRUE, so this means this file contains only motif hits by variants

# how many variants hits contain motifs in each allele
table((variant_hits %>% group_by(variant_loc, motif_name) %>% summarise(n_alleles=n_distinct(allele)))$n_alleles)
# ref 15647  alt 2623 

#get motif info
### modisco patterns
modisco_out <- read_tsv("~/cd4_chrombpnet/chrombpnet_model_b7/tfmodisco_motifs_hocomoco_jaspar_cisbp/model/motifs_with_tf_annotated.tsv") 

modisco_out<- modisco_out %>%      # drop any stray spaces / tabs
  separate(pattern,
           into = c("class", "pattern_id"),
           sep  = "\\.",                     # split at the dot
           remove = FALSE) %>% mutate(TF_best_match_code_friendly_unique=make.unique(modisco_out$TF_best_match_code_friendly)) %>%
  mutate(match_confidence = case_when(
    str_starts(TF_best_match_code_friendly_unique, "UNKNOWN")  ~ "no_TF_match",
    str_starts(TF_best_match_code_friendly_unique, "PUTATIVE") ~ "weak_TF_match",
    TRUE                              ~ "strong"
  )) %>% group_by(TF_best_match_code_friendly_unique) %>%
  mutate(dup_count = row_number()) %>%
  ungroup() %>%
  mutate(
    TF_best_match_code_friendly_unique_unique = if_else(
      duplicated(TF_best_match_code_friendly_unique) | duplicated(TF_best_match_code_friendly_unique, fromLast = TRUE),
      paste0(TF_best_match_code_friendly_unique, "_", dup_count),
      TF_best_match_code_friendly_unique
    )
  ) %>%
  dplyr::select(-dup_count)
head(modisco_out)

## annotate motif family
# --- 1) Canonical family map (named character vector) ---
tf_map <- c(
  "ATF1_CREB1"              = "AP-1/ATF/CREB (bZIP)",
  "CTCF"                    = "CTCF",
  "PAX5"                    = "PAX",
  "KLF12"                   = "KLF (C2H2 ZF)",
  "RUNX1_BCL11B"            = "Mixed: RUNX/BTB-ZF",
  "ETS1"                    = "ETS",
  "TYY1"                    = "YY (C2H2 ZF)",  
  "NRF1"                    = "NRF (bZIP)",
  "NFYA_B_C"                = "NF-Y (CCAAT)",
  "PUTATIVE_RUNX2"          = "RUNX",
  "JUN_ATF"                 = "AP-1/ATF/CREB (bZIP)",
  "UNKNOWN_2_CTCF_LIKE"     = "CTCF",
  "IRF1_IRF2"               = "IRF",
  "ZBTB33_KAISO"            = "BTB/POZ ZF",
  "CREB1_AP1"               = "AP-1/ATF/CREB (bZIP)",
  "PUTATIVE_ETS1_IKZF3_ERG" = "Mixed: ETS/IKAROS",
  "ZNF76_ZN143_1"           = "ZNF (C2H2 ZF)",
  "PUTATIVE_RUNX1"          = "RUNX",
  "TFE3"                    = "MiT/TFE (bHLH-ZIP)",
  "PUTATIVE_ETS_IRF"        = "Mixed: ETS/IRF",
  "ELF5"                    = "ETS",
  "IRF4_IRF8"               = "IRF",
  "RFX3_RFX5_RFX1"          = "RFX",
  "PUTATIVE_RORG"           = "Nuclear receptor (ROR)",
  "TCF7_LEF1"               = "TCF/LEF (HMG)",
  "IRF9"                    = "IRF",
  "CTCF.1"                  = "CTCF",
  "ZNF76_ZN143_2"           = "ZNF (C2H2 ZF)",
  "NFKB1"                   = "NF-κB (Rel)",
  "RFX5"                    = "RFX",
  "IRF1"                    = "IRF",
  "BCL11B"                  = "BTB/POZ ZF",
  "UNKNOWN_1_ETS_LIKE"      = "ETS",
  "SP2_SP3"                 = "Sp/KLF (C2H2 ZF)",
  "KLF6"                    = "KLF (C2H2 ZF)"
)

# Join the updated naming annotation to the hits
tf_map_tbl <- enframe(tf_map, name = "motif", value = "family")

motif_col <- "TF_best_match_code_friendly_unique_unique"   # <-- column in `unique_hits` that has the motif names

modisco_out <- modisco_out %>%
  mutate(.motif_key = str_trim(.data[[motif_col]])) %>%        # trim whitespace just in case
  left_join(tf_map_tbl, by = c(".motif_key" = "motif")) %>%    # add `family`
  dplyr::select(-.motif_key)

head(modisco_out)


################
###############
## The filtering motif motif patterns and their hits needs to follow a few rules:
## - (1) exclude motif patterns that do not have at least 40 seqlets to back them up
## - (1) de-duplication of overlapping patterns instances per motif with a 3bp overlap
## - (3) keep only hits with a hit_correlation > 0.8


#-- annotate motf variant instances with family name
motif_ids<- modisco_out %>% dplyr::select(pattern, num_seqlets, match_confidence, TF_best_match_code_friendly_unique_unique, family)
head(variant_hits)
variant_hits<-variant_hits %>% left_join(motif_ids, join_by("motif_name"=="pattern"))


## (1) filter per min seqlets
modisco_out40 <- modisco_out %>% filter(num_seqlets>40) #motif version
write_tsv(modisco_out40, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/modisco_motids_minseq40.tsv")

print(paste0("After filtering for motifs patterns with higher than 40 seqlets to back them up, there are ", n_distinct(modisco_out40$family), " remaining Motif families"))

variant_hits40 <- variant_hits %>% filter(num_seqlets>40) #instance version

## qc

metrics  <- c("hit_coefficient", "hit_correlation", "hit_importance")
filters  <- c(5, 0.8, 0.1)   # same order as metrics

pt.list <- vector("list", length(metrics))
names(pt.list) <- metrics

for (i in seq_along(metrics)) {
  metric  <- metrics[i]
  cutoff  <- filters[i]
  
  label_text <- paste0(metric, " = ", cutoff)
  
  p <- ggplot(variant_hits40, aes(x = .data[[metric]])) +
    geom_histogram(bins = 30, fill = "steelblue2", color = "gray20", linewidth = 0.2) +
    geom_vline(xintercept = cutoff,
               linetype  = "dashed",
               linewidth = 0.6,
               color     = "black") +
    annotate(
      "text",
      x     = cutoff,
      y     = Inf,              
      label = label_text,
      vjust = 2,              
      size  = 4
    ) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = 16),
      axis.title = element_text(size = 16)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(metric) +
    ylab("Number of Instances")
  
  pt.list[[i]] <- p
}

pt.filters<-plot_grid(plotlist = pt.list, ncol = 3)
ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/variant_motif_overlap/motif_instance_qc.pdf", width = 12, height = 3,  dpi=300)


## (1) filter seqets
variant_hits40_coeff <- variant_hits %>% filter(hit_coefficient>5 & hit_correlation>0.8 & hit_importance>0.1) #instance version

    # pt.list_filt <- vector("list", length(metrics))
    # names(pt.list_filt) <- metrics
    # 
    # for (i in seq_along(metrics)) {
    #   metric <- metrics[i]
    #   
    #   p <- ggplot(variant_hits40_coeff, aes(x = .data[[metric]])) +
    #     geom_histogram(bins = 30, fill = "steelblue2", color = "gray20", linewidth = 0.2) +
    #     theme_classic() +
    #     theme(
    #       axis.text  = element_text(size = 16),
    #       axis.title = element_text(size = 16)
    #     ) +
    #     scale_y_continuous(expand = c(0, 0)) +
    #     xlab(metric) +
    #     ylab("Number of Instances")
    #   
    #   pt.list_filt[[i]] <- p
    # }
    # 
    # plot_grid(plotlist = pt.list_filt, ncol = 3)  

#number of motif hits per allele
table(variant_hits40_coeff$allele)

#how many distinct variants overlap a motif?
n_distinct(variant_hits40_coeff$variant_loc)

# how many motifs are present in both alleles
table((variant_hits40_coeff %>% group_by(variant_loc, motif_name) %>% summarise(n_alleles=n_distinct(allele)))$n_alleles)

variant_summary <- variant_hits40_coeff %>%
    group_by(peak_id, variant_loc, family) %>%
  summarise(n_alleles = n_distinct(allele), .groups = "drop") %>%
  mutate(
    Disruption_level = ifelse(n_alleles == 2, "mild", "strong")
  )

family_counts <- variant_hits40_coeff %>%
  group_by(variant_loc, allele, family) %>%
  summarise(
    n_hits          = n(),                              # total rows (hits)
    n_unique_motifs = n_distinct(motif_name),           # distinct motifs
    motifs_collapsed = paste(sort(unique(TF_best_match_code_friendly_unique_unique)),  # collapse names
                             collapse = ";"),
    .groups = "drop"
  )

# 2) Any variant_loc × allele × family with >1 *distinct* motif_name?
family_counts_multi_motifs <- family_counts %>%
  filter(n_unique_motifs > 1)

## (2) deduplicate instances

# Overlapping hits (≥4 bp overlap) are clustered together across all motifs in that family.
# Each cluster is represented by a single row: the hit with max hit_correlation.
# motifs_collapsed on that row gives you all motif names in that overlapping cluster, e.g.
# IRF1;IRF2;IRF4.

# convert to granges
gr <-GRanges(
    seqnames = variant_hits40_coeff$chr,
    ranges   = IRanges(start = variant_hits40_coeff$start, end = variant_hits40_coeff$end),
    motif_name      = variant_hits40_coeff$TF_best_match_code_friendly_unique_unique,
    pattern_name    = variant_hits40_coeff$motif_name,
    family          = variant_hits40_coeff$family,
    hit_correlation = variant_hits40_coeff$hit_correlation,
    variant_loc = variant_hits40_coeff$variant_loc,
    allele = variant_hits40_coeff$allele,
    hit_coefficient = variant_hits40_coeff$hit_coefficient,
    match_confidence= variant_hits40_coeff$match_confidence,
    hit_importance = variant_hits40_coeff$hit_importance,
    num_seqlets = variant_hits40_coeff$num_seqlets,
    peak_id=variant_hits40_coeff$peak_id
  )


#------------------------------------------------------------
# helper: keep top-correlated hit per ≥min_olap cluster,
#         and record all motif_names in that cluster
filter_per_family <- function(gr_family, min_olap = 4L) {
  n <- length(gr_family)
  if (n == 0L) return(gr_family)
  
  # Start with each hit having its own motif_name as collapsed label
  mcols(gr_family)$motifs_collapsed <- as.character(mcols(gr_family)$motif_name)
  
  # Overlaps (≥ min_olap bp)
  ov <- findOverlaps(gr_family, minoverlap = min_olap, ignore.strand = TRUE)
  
  # If no overlaps at all in this family, nothing to dedup
  if (length(ov) == 0L) {
    return(gr_family)
  }
  
  # Build graph: nodes = hits (by index), edges = overlaps
  g <- graph_from_edgelist(
    cbind(queryHits(ov), subjectHits(ov)),
    directed = FALSE
  )
  
  comps <- components(g)$membership        # component ID per vertex
  vertex_ids <- as.integer(names(comps))   # these are the original row indices in gr_family
  
  keep_idx <- integer(0)
  
  # For each connected component (cluster of overlapping hits):
  for (comp_id in unique(comps)) {
    idx <- vertex_ids[comps == comp_id]   # row indices in this cluster
    
    # pick the winner = max hit_correlation
    best <- idx[which.max(mcols(gr_family)$hit_correlation[idx])]
    
    # collapsed motif names within this cluster
    collapsed <- paste(
      sort(unique(mcols(gr_family)$motif_name[idx])),
      collapse = ";"
    )
    
    mcols(gr_family)$motifs_collapsed[best] <- collapsed
    
    keep_idx <- c(keep_idx, best)
  }
  
  # Hits that *never* overlapped anyone (e.g. short hits < min_olap, etc.)
  not_in_graph <- setdiff(seq_len(n), vertex_ids)
  if (length(not_in_graph) > 0L) {
    # They already have motifs_collapsed = motif_name, so just keep them
    keep_idx <- c(keep_idx, not_in_graph)
  }
  
  # Return only the chosen rows, sorted by genomic order within this family
  gr_family[sort(keep_idx)]
}

#-deduplicate motifs
filtered_gr_family <- gr %>%
  split(list(
    variant_loc = mcols(.)$variant_loc,
    family      = mcols(.)$family,
    allele      = mcols(.)$allele
  ), drop = TRUE) %>%
  lapply(filter_per_family, min_olap = 4L) %>%
  GRangesList() %>%
  unlist(use.names = FALSE)


filtered_gr_family

#how many distinct variants overlap a motif?
n_distinct(filtered_gr_family$variant_loc)
n_distinct(filtered_gr_family$family)



saveRDS(filtered_gr_family, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/filtered_motif_instances.rds")

#save df version
df_motifs <- as.data.frame(filtered_gr_family)
write_tsv(df_motifs, "~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/filtered_variant_hits.tsv")

filtered_gr_family<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/filtered_variant_hits.tsv")
unique_hits<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap_hit_caller/unique_hits_filtered_motif_instances.tsv")


filtered_gr_family[filtered_gr_family$ranges==150567551-150567559]
