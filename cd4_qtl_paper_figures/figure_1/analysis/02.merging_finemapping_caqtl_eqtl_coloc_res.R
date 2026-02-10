.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(GenomicRanges)
library(cowplot)
library(LaCroixColoR)
library(rcartocolor)
library(pheatmap)
library(stringr)
options(bitmapType = "cairo")
source("/gchm/cd4_qtl_paper_figures/figure_1/helper_functions.R")
library(BiocParallel)
register(MulticoreParam(workers = 8))  # or SnowParam() for SLURM

# File paths
coloc_files <- system("ls /gcgl/sghatan/marlis_pj/coloc/coloc_results/ca_eqtl_coloc/*_coloc_results.csv", intern = TRUE)
eqtl_finemapping_files <- system("ls /gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_credible_sets/All_CD4T_cells/All_CD4T_cells_chr*_credible_sets.txt", intern = TRUE)
caqtl_finemapping_files <- system("ls /gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_credible_sets/CD4T_chromatin/CD4T_chromatin_chr*_credible_sets.txt", intern = TRUE)
peak_count_path <- "/gchm/ATAC-seq_analysis/diff_accesibility_ana/results/peak_counts/RAW_cd4_atac_peakcounts_ranges_scrna_union.csv"
correlation_results <- read.csv("~/cd4_qtl_paper_figures/figure_1/data/caPeak_coacc_correlation_results_allchrs.csv")

# Load peak data
ranges.table <- read.csv2(peak_count_path, sep = ",") %>%
  filter(Chr %in% paste0("chr", 1:22))
rownames(ranges.table) <- ranges.table$X

gr <- GRanges(seqnames = ranges.table$Chr,
              ranges = IRanges(ranges.table$Start, ranges.table$End),
              mcols = data.frame(peakID = rownames(ranges.table)))
names(gr) <- gr$mcols.peakID

# Load correlated peak data
sig_corrs <- correlation_results %>% filter(FDR < 0.05, cor > 0)
lead_to_corr <- split(sig_corrs$corPeak, sig_corrs$leadPeak)
lead_corr_df <- sig_corrs %>% select(leadPeak, corPeak)
target_gr <- gr[names(gr) %in% unique(lead_corr_df$corPeak)]


# Extract chromosome labels
chr <- str_extract(coloc_files, "[^_/]+(?=_coloc_results)")

# Define summarization function
summarise_coloc <- function(x, gr, lead_corr_df, lead_to_corr, target_gr) {
  message("Processing chromosome: ", chr[x])
  
  coloc <- suppressMessages(read_delim(coloc_files[x]))
  eqtl_finemapping <- suppressMessages(read_delim(eqtl_finemapping_files[x])) %>% prep_finemap_df()
  caqtl_finemapping <- suppressMessages(read_delim(caqtl_finemapping_files[x])) %>% prep_finemap_df()
  
  ### Classify caQTL variants into C1–C4
  # variant-level classification 
  # ─────────────────────────────────────────────────────────────
  # Variant-level peak classification (Category):
  # Each fine-mapped variant is assigned a category based on its 
  # overlap with chromatin accessibility peaks:
  #   - C1: Variant overlaps its assigned peak
  #   - C2: Variant does not overlap assigned peak, but overlaps a peak 
  #         correlated with the assigned peak
  #   - C3: Variant does not overlap assigned or correlated peaks, 
  #         but overlaps another open peak
  #   - C4: Variant does not overlap any peak
  # This is used to infer the regulatory relevance of variants 
  # relative to their local chromatin context.
  # ─────────────────────────────────────────────────────────────
  caqtl_finemapping$variant_pos <- as.numeric(sub(".*:(\\d+)\\[.*", "\\1", caqtl_finemapping$variant_id))
  
  variantsGR <- GRanges(seqnames = caqtl_finemapping$chr,
                        ranges = IRanges(start = caqtl_finemapping$variant_pos, end = caqtl_finemapping$variant_pos),
                        assigned_peak = caqtl_finemapping$peak)
  
  peak_names <- unique(caqtl_finemapping$peak)
  caQTL_peaksGR <- gr[names(gr) %in% peak_names]
  overlaps <- findOverlaps(variantsGR, caQTL_peaksGR)
  
  caqtl_finemapping$InAssignedPeak <- FALSE
  caqtl_finemapping$InNearbyPeak <- FALSE
  
  for (i in seq_along(queryHits(overlaps))) {
    snp_idx <- queryHits(overlaps)[i]
    overlapping_peak <- names(caQTL_peaksGR)[subjectHits(overlaps)[i]]
    assigned_peak <- caqtl_finemapping$peak[snp_idx]
    if (overlapping_peak == assigned_peak) {
      caqtl_finemapping$InAssignedPeak[snp_idx] <- TRUE
    } else {
      caqtl_finemapping$InNearbyPeak[snp_idx] <- TRUE
    }
  }
  
  caqtl_finemapping$InCorrelatedPeak <- FALSE
  snp_idx <- which(caqtl_finemapping$peak %in% lead_corr_df$leadPeak)
  
  for (i in snp_idx) {
    lead <- caqtl_finemapping$peak[i]
    corr_peaks <- lead_to_corr[[lead]]
    target_subset <- target_gr[names(target_gr) %in% corr_peaks]
    if (length(target_subset) == 0) next
    single_snp <- variantsGR[i]
    if (length(findOverlaps(single_snp, target_subset)) > 0) {
      caqtl_finemapping$InCorrelatedPeak[i] <- TRUE
    }
  }
  
  caqtl_finemapping$Category <- "C4"
  caqtl_finemapping$Category[caqtl_finemapping$InAssignedPeak] <- "C1"
  caqtl_finemapping$Category[!caqtl_finemapping$InAssignedPeak & caqtl_finemapping$InCorrelatedPeak] <- "C2"
  caqtl_finemapping$Category[!caqtl_finemapping$InAssignedPeak & !caqtl_finemapping$InCorrelatedPeak & caqtl_finemapping$InNearbyPeak] <- "C3"
  
  # CS-level (credible set) categorization
  # ─────────────────────────────────────────────────────────────
  # CS-level categorization (cs_type_any_var):
  # Each fine-mapped credible set (CS) is assigned a category based on
  # the hierarchical presence of its variants relative to chromatin peaks:
  #   - C1 ("in_caPeak"): ≥1 variant overlaps its assigned peak
  #   - C2 ("in_corr_Peak"): no variant in assigned peak, but ≥1 in correlated peak
  #   - C3 ("in_other_Peak"): not C1/C2, but ≥1 in any other peak
  #   - C4 ("no_Peak_overlap"): no variant in any peak
  # The cs_type_any_var field reflects this hierarchy at the CS level.
  # ─────────────────────────────────────────────────────────────
  
  cs_category_presence <- caqtl_finemapping %>%
    group_by(finemapped_cs, Category) %>%
    summarise(has = any(TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Category, values_from = has, values_fill = FALSE)
  
  cs_category_presence <- cs_category_presence %>%
    mutate(cs_type_any_var = case_when(
      C1                           ~ "in_caPeak",
      !C1 & C2                     ~ "in_corr_Peak",
      !C1 & !C2 & C3               ~ "in_other_Peak",
      TRUE                         ~ "no_Peak_overlap"
    ))
  
  caqtl_finemapping <- left_join(caqtl_finemapping, cs_category_presence, by = "finemapped_cs")
  
  
  ### Merge with coloc results
  coloc <- coloc %>%
    mutate(
      finemapped_cs_caqtl = paste0(peak, "_", region.x, "_", idx1),
      finemapped_cs_eqtl  = paste0(gene, "_", region.y, "_", idx2)
    )
  
  coloc_simple <- coloc %>% 
    select(caQTL_variant, eQTL_variant, starts_with("PP."), 
           finemapped_cs_caqtl, finemapped_cs_eqtl)
  
  caqtl_coloc <- left_join(caqtl_finemapping, coloc_simple, 
                           join_by(finemapped_cs == finemapped_cs_caqtl)) %>%
    mutate(finemapped_cs_caqtl = finemapped_cs)
  
  allres <- full_join(caqtl_coloc, eqtl_finemapping, 
                      join_by(finemapped_cs_eqtl == finemapped_cs)) %>%
    mutate(
      coloc_status = case_when(
        !is.na(finemapped_cs_caqtl) & !is.na(finemapped_cs_eqtl) ~ "coloc",
        !is.na(finemapped_cs_caqtl) & is.na(finemapped_cs_eqtl)  ~ "caQTL_only",
        is.na(finemapped_cs_caqtl) & !is.na(finemapped_cs_eqtl)  ~ "eQTL_only",
        TRUE                                                     ~ "unassigned"
      ),
      finemapped_cs_coloc = case_when(
        !is.na(finemapped_cs_caqtl) & !is.na(finemapped_cs_eqtl) ~ paste0(finemapped_cs_caqtl, "_", finemapped_cs_eqtl),
        !is.na(finemapped_cs_caqtl) & is.na(finemapped_cs_eqtl)  ~ finemapped_cs_caqtl,
        is.na(finemapped_cs_caqtl) & !is.na(finemapped_cs_eqtl)  ~ finemapped_cs_eqtl,
        TRUE                                                     ~ NA_character_
      ),
      chromosome = chr[x]
    )
  # Write to file
  out.dir<-"/gcgl/mmatos/cd4_aging_project/cd4_qtl_data/colocalization/caqtl_eqtl_coloc_finemapped/"
  allres_outfile <- file.path(out_dir, paste0("eqtl_caqtl_finemap_coloc_summary_chr", chr[x], ".tsv.gz"))

  write_tsv(allres, allres_outfile)
  return(list(allres = allres, coloc_simple = coloc_simple))
}


# Run over chromosomes
results <- bplapply(seq_along(coloc_files), function(x) {
  summarise_coloc(x, gr, lead_corr_df, lead_to_corr, target_gr)
})


# Combine all outputs
summary_df       <- bind_rows(lapply(results, `[[`, "allres"))
coloc_simple_all <- bind_rows(lapply(results, `[[`, "coloc_simple"))


#a credible-set level ccolocalization summary between peaks and genes
write_delim(summary_df, "/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv", delim = "\t") # a variant level colocalization (all finemappped)
write_delim(coloc_simple_all, "/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_coloc_allchrs.tsv", delim = "\t") #a credible-set level ccolocalization summary between peaks and genes


print("Finished!")

# out.dir<-"/gcgl/mmatos/cd4_aging_project/cd4_qtl_data/colocalization/caqtl_eqtl_coloc_finemapped/"
# 
# summary_df %>%
#   group_by(chromosome) %>%
#   group_walk(~ write_tsv(.x, file = file.path(out.dir, paste0("eqtl_caqtl_finemap_coloc_summary_", .y$chromosome, ".tsv"))))


