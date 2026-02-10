# -------------------------------
# Mirror plot: Forward (peak->gene) above, Reverse (gene->peak) below
# Requires: forward/reverse TSVs already exist (concatenated)
# -------------------------------

.libPaths(c(
  "/gchm/R/x86_64-pc-linux-gnu-library/4.4",
  "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"
))
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  library(rtracklayer)
  library(purrr)
  library(igraph)
  library(readr)
  library(cowplot)
})

# -------------------------------
# Inputs
# -------------------------------
forward_fn <- "~/cd4_qtl_paper_figures/figure_3/data/mediation_res_filter_cis_only/concatenated/caqtl_mediation_cis_annotated_allchr.tsv"
reverse_fn <- "~/cd4_qtl_paper_figures/figure_3/data/mediation_res_filter_cis_only/concatenated/eqtl_mediation_cis_annotated_allchr.tsv"

chrom_sizes_path <- "~/resources/genome/hg38.p14.chrom.sizes"
peak_coords_path <- "~/cd4_caQTL_analysis/variant_to_peak_QTL/run_012625_qc_aware_qsmooth_CPM_MAF5_FDR5_1MB/results/003_inputs/filtered_qsmooth_norm_cpm/cd4_atac_processed_peaks_coordinates.bed"
gtf_path <- "resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.basic.annotation.gtf.gz"

# Mirror plot parameters
baseline <- 0.5
cutoff   <- NULL
n_labels_forward <- 40  # gene labels (top)
n_labels_reverse <- 20  # peak labels (bottom)

# -------------------------------
# Load data
# -------------------------------
forward <- fread(forward_fn)
reverse <- fread(reverse_fn)

# -------------------------------
# Standardize dfs to merge
# -------------------------------

columnnames<-c("finemapped_cs_coloc", "exposure_variable" ,  "chromosome"  ,  "phenotype_A"    ,   "phenotype_B"    ,   "coloc_triplet" ,  "phenotype_B_instrument_proximity" ,  "p2" , "p3" ,  "p4",    
               "p5",   "p2p3"  ,  "fdr_sig_1", "p2p3_cutoff_fdr_1"  , "fdr_sig_5"  , "p2p3_cutoff_fdr_5")
colnames(forward) <-columnnames
colnames(reverse) <-columnnames

# Common columns, ordered like iv_5rev
common_cols <- names(reverse)[names(reverse) %in% names(forward)]

# Subset & reorder both data frames
reverse <- reverse %>%
  dplyr::select(all_of(common_cols)) %>%
  mutate(type = "reverse",
         variant_position = as.integer(sub(".*:(\\d+)\\[.*", "\\1", exposure_variable))) # Extract variant position from "1:16504530[b38]CTT,CT"

forward <- forward %>%
  dplyr::select(all_of(common_cols)) %>%
  mutate(type = "forward",
         variant_position = as.integer(sub(".*:(\\d+)\\[.*", "\\1", exposure_variable)))

# -------------------------------
# concatenate findr results for both models
# -------------------------------

findr_res<-bind_rows(forward, reverse) #these are all the results in cis, not filtered for significance
findr_res<-findr_res %>% mutate(causality_uniq=paste0( type, "_", finemapped_cs_coloc, "_", exposure_variable, "_", phenotype_A, "_", phenotype_B))

# -------------------------------
# collapse findr res by coloc locus
#
# Since coloc loci are interrelated, there is a loci of overlap per loci
# i am doing a three way association, meaning that if A -> B and B->C then A->C
# -------------------------------

#1 Take one row per unique triad
unique_triads <- findr_res %>%
  distinct(finemapped_cs_coloc, .keep_all = TRUE)
#2 Extract the two CS IDs per triad
#cs_pattern <- "[0-9XY]+:[0-9]+-[0-9]+_[0-9]+"
cs_pattern <- "(?<=_|^)[A-Za-z0-9_.-]+_[0-9XY]+:[0-9]+-[0-9]+_[0-9]+"

triad_cs <- unique_triads %>%
  mutate(
    cs_ids = str_extract_all(finemapped_cs_coloc, cs_pattern)
  )

# Check how many matches per row
table(lengths(triad_cs$cs_ids))

# sanity: each triad should have exactly 2 CS IDs
stopifnot(all(lengths(triad_cs$cs_ids) == 2))

triad_cs <- triad_cs %>%
  mutate(
    cs1 = map_chr(cs_ids, 1),
    cs2 = map_chr(cs_ids, 2)
  )

#3 Build the CS–CS graph and get locus IDs
# Edges between CSs from each unique triad
edges <- triad_cs %>%
  dplyr::select(cs1, cs2) %>%
  distinct()

g <- graph_from_data_frame(edges, directed = FALSE)

comp <- components(g)

rm(g)
rm(edges)

cs_to_locus <- tibble(
  cs_id    = names(comp$membership),
  locus_id = paste0("locus_", comp$membership)
)

rm(comp)
#4 Assign each triad a locus_id
triad_to_locus <- triad_cs %>%
  left_join(cs_to_locus, by = c("cs1" = "cs_id")) %>%
  dplyr::select(finemapped_cs_coloc, cs1, cs2, locus_id)

#5 Map loci back to the full findr_res (with duplicates)
findr_res <- findr_res %>%
  left_join(
    triad_to_locus %>% select(finemapped_cs_coloc, locus_id),
    by = "finemapped_cs_coloc"
  )
rm(triad_to_locus)
rm(forward)
rm(reverse)

write_tsv(findr_res, "~/cd4_qtl_paper_figures/figure_3/data/findr_MEDIATION_All_results_collapsed_not_filtered.tsv")

# # caPeak name = first 5 underscore tokens
# forward <- forward %>%
#   mutate(
#     caPeak_name = str_replace(
#       finemapped_cs_coloc,
#       "^(([^_]+_){4}[^_]+).*",
#       "\\1"
#     )
#   )

# -------------------------------
# Chrom sizes -> chr offsets
# -------------------------------
chrom_sizes_dt <- fread(chrom_sizes_path, col.names = c("chr", "chr_len"))

chr_tbl <- as.data.frame(chrom_sizes_dt) %>%
  transmute(chr = chr, chr_len = as.numeric(chr_len)) %>%
  filter(grepl("^chr([1-9]|1[0-9]|2[0-2])$", chr)) %>%  # chr1..chr22 only
  mutate(chromosome = as.integer(sub("^chr", "", chr))) %>%
  arrange(chromosome) %>%
  mutate(
    chr_start = c(0, cumsum(chr_len)[-n()]),
    chr_mid   = chr_start + chr_len / 2
  ) %>%
  select(chromosome, chr_len, chr_start, chr_mid)

# -------------------------------
# Peak coordinates + summit
# -------------------------------
peak_coords <- fread(peak_coords_path)
colnames(peak_coords) <- c("peak_chr", "peak_start", "peak_end", "peak_name")

forward <- findr_res %>% filter(type=="forward") %>%
  left_join(peak_coords, by = join_by(phenotype_A == peak_name)) %>%
  mutate(
    peak_summit = as.integer(round(peak_start + (peak_end - peak_start) / 2))
  )

# Peak genome position for forward plot
forward_plot <- forward %>%
  mutate(chromosome = as.integer(sub("^chr", "", chromosome))) %>%
  left_join(chr_tbl %>% select(chromosome, chr_start), by = "chromosome") %>%
  mutate(peak_genome_pos = chr_start + peak_summit) %>%
  select(locus_id,peak_genome_pos, p2p3, fdr_sig_5, coloc_triplet, phenotype_A, phenotype_B) %>%
  distinct()

# -------------------------------
# Gene TSS from GTF (primary assembly)
# -------------------------------
gtf <- import(gtf_path)

gene_tss <- as.data.frame(gtf[gtf$type == "gene"]) %>%
  transmute(
    gene_id   = sub("\\..*$", "", gene_id),
    gene_name = gene_name,
    chr       = as.character(seqnames),
    strand    = as.character(strand),
    tss       = if_else(strand == "+", start, end)
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%      # drop alt contigs
  arrange(gene_name, chr, tss) %>%
  group_by(gene_name) %>%
  slice_head(n = 1) %>%
  ungroup()

# -------------------------------
# Reverse: join TSS + chr offsets + (optional) peak coords for labels
# -------------------------------
reverse <- findr_res %>% filter(type=="reverse") %>%
  mutate(phenotype_A = sub("\\..*$", "", phenotype_A)) %>%
  left_join(gene_tss, join_by("phenotype_A"=="gene_name")) %>%
  filter(!is.na(tss), !is.na(chr)) %>%
  mutate(chromosome = as.integer(sub("^chr", "", chr))) %>%
  left_join(chr_tbl %>% select(chromosome, chr_start), by = "chromosome") %>%
  mutate(gene_genome_pos = chr_start + tss) %>%
  left_join(peak_coords, join_by("phenotype_B"=="peak_name")) %>%
  mutate(
    peak_coord = paste0(peak_chr, ":", peak_start, "-", peak_end)
  )

reverse_plot <- reverse %>%
  select(locus_id,gene_genome_pos, p2p3, fdr_sig_5, coloc_triplet, phenotype_A, phenotype_B, peak_coord) %>%
  distinct()

# -------------------------------
# Build mirrored dataset
#   forward: y = +(p2p3 - baseline)
#   reverse: y = -(p2p3 - baseline)
# -------------------------------
mir_u <- mir %>%
  distinct(dir, genome_pos, y, label, .keep_all = TRUE)

mir_f <- forward_plot %>%
  transmute(
    dir = "forward",
    genome_pos = peak_genome_pos,
    y =  p2p3,
    p2p3,
    fdr_sig_5,
    coloc_triplet,
    label = phenotype_B
  )

mir_r <- reverse_plot %>%
  transmute(
    dir = "reverse",
    genome_pos = gene_genome_pos,
    y = -p2p3,
    p2p3,
    fdr_sig_5,
    coloc_triplet,
    label = peak_coord   # or peak_name
  )

mir <- bind_rows(mir_f, mir_r)
ymax <- max(abs(mir$y), na.rm = TRUE)   # should be <= 1

# Label top hits per side (one label per unique label)
labs_f <- mir_u %>%
  filter(dir == "forward") %>%
  group_by(label) %>%
  slice_max(p2p3, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(p2p3, n = n_labels_forward, with_ties = FALSE)

labs_r <- mir_u %>%
  filter(dir == "reverse") %>%
  group_by(label) %>%
  slice_max(p2p3, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  slice_max(p2p3, n = n_labels_reverse, with_ties = FALSE)

labs_mir <- mir %>%
  group_by(dir, label) %>%
  slice_max(p2p3, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(dir) %>%
  arrange(dplyr::desc(p2p3), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  filter((dir == "forward" & rank <= n_labels_forward) |
           (dir == "reverse" & rank <= n_labels_reverse)) %>%
  select(-rank)

library(scales)

k <- 0.1  # <--- smaller = more compression of 0–0.5

compress_0_5 <- trans_new(
  name = "compress_0_5",
  transform = function(y) {
    s <- sign(y); a <- abs(y)
    a_t <- ifelse(a <= 0.5, a * k, 0.5 * k + (a - 0.5))
    s * a_t
  },
  inverse = function(yt) {
    s <- sign(yt); a <- abs(yt)
    a_inv <- ifelse(a <= 0.5 * k, a / k, 0.5 + (a - 0.5 * k))
    s * a_inv
  }
)

# -------------------------------
# Genes to highlight
# -------------------------------
#unique((forward_plot %>% filter(fdr_sig_5==TRUE & coloc_triplet==TRUE))$phenotype_B)
tcell_like <- c(
  "TRBC1","TRBV4-2","TRAV8-5","KIR3DL2","ULBP3",
  "TNFRSF18","IL23R","PTGER4","CD84","CD300A","LAIR2",
  "BACH2","HIVEP3","DOCK8","DOCK7","HOPX","GNAI1","ARL4A",
  "ARHGAP5","PLEKHG1","CHN1","FAM13A",
  "FAS","IL32","ISG15","MX1","IFI27L2","LY6E-DT",
  "GSDMD","NOD2","TLR1","APOBEC3G","SECTM1",
  "TNFSF13B","ALOX5AP","ORMDL3","GCNT1","VAMP5","CMTM8","LY86","LY86-AS1"
)
labs_tcell <- mir %>%
  filter(dir == "forward", label %in% tcell_like) %>%
  group_by(label) %>%
  slice_max(p2p3, n = 1, with_ties = FALSE) %>%
  ungroup()
labs_mir2 <- labs_mir %>%
  filter(dir != "forward")   # keep only reverse labels from the generic set
mir <- mir %>%
  mutate(fill_dir_sig = ifelse(fdr_sig_5, as.character(dir), "not_sig"))
# -------------------------------
# Plot
# -------------------------------
p <- ggplot(mir, aes(
  x = genome_pos, y = y,
  fill  = fill_dir_sig,
  shape = factor(coloc_triplet),
  size  = fdr_sig_5
)) +
  geom_vline(
    data = chr_tbl, aes(xintercept = chr_start),
    inherit.aes = FALSE, linewidth = 0.2, alpha = 0.3
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  
  # points
  geom_point(alpha = 0.7, stroke = 0.1, color = "grey25") +
  
  # reverse labels (peaks) — explicitly anchor to their x/y
  geom_label_repel(
    data = labs_mir2,
    inherit.aes = FALSE,
    aes(x = genome_pos, y = y, label = label),
    size = 3, nudge_y = -0.06,
    box.padding = 0.3, point.padding = 0.2,
    label.size = 0.25, fill = "white", alpha = 0.95,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  
  # T-cell gene labels (forward) — explicitly anchor to their x/y
  geom_label_repel(
    data = labs_tcell,
    inherit.aes = FALSE,
    aes(x = genome_pos, y = y, label = label),
    size = 3.2, nudge_y = 0.08,
    box.padding = 0.35, point.padding = 0.25,
    label.size = 0.25, fill = "white", alpha = 0.98,
    max.overlaps = Inf, min.segment.length = 0
  ) +
  
  scale_x_continuous(
    breaks = chr_tbl$chr_mid,
    labels = chr_tbl$chromosome,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans  = compress_0_5,
    limits = c(-1.1, 1.1),
    breaks = c(-1, -0.9, -0.5, 0, 0.5, 0.9, 1),
    labels = function(v) sprintf("%.1f", abs(v)),
    expand = c(0, 0)
  ) +
  
  scale_fill_manual(
    values = c(not_sig = "grey90", forward = "#9ccb86", reverse = "#ef7fc5"),
    drop = FALSE,
    name = "Significant (FDR 5%)"
  ) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), name = "Coloc triplet") +
  scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.5), guide = "none") +
  labs(x = "Chromosome", y = "p2p3 (back-to-back)") +
  theme_classic()

print(p)


# Optional save
ggsave("~/cd4_qtl_paper_figures/figure_3/plots/dec2025/mediation_posterior_mirror_plot.pdf", p, width = 12, height = 5, useDingbats = FALSE)

##plots
findr_res<-fread("~/cd4_qtl_paper_figures/figure_3/data/findr_MEDIATION_All_results_collapsed_not_filtered.tsv")
totals<-fread("~/cd4_qtl_paper_figures/figure_3/data/findr_MEDIATION_inputs_collapsed.tsv") 
findr_res_sig<-findr_res %>% filter(fdr_sig_5==TRUE)

res.base <- findr_res_sig %>%
  mutate(
    triads = paste(locus_id, phenotype_A, phenotype_B, sep = "_")
  ) %>%
  distinct(type, triads, locus_id, phenotype_A, phenotype_B, fdr_sig_5, fdr_sig_1)

res_5 <- res.base %>%
  filter(fdr_sig_5 == TRUE) %>%          # all 5% hits, including 1%
  group_by(type) %>%
  summarise(
    n_triads    = n_distinct(triads),
    n_loci      = n_distinct(locus_id),
    n_mediators = n_distinct(phenotype_A),
    n_targets   = n_distinct(phenotype_B),
    .groups     = "drop"
  ) %>%
  mutate(significance = "5% FDR")
res_1 <- res.base %>%
  filter(fdr_sig_1 == TRUE) %>%         # only 1% hits
  group_by(type) %>%
  summarise(
    n_triads    = n_distinct(triads),
    n_loci      = n_distinct(locus_id),
    n_mediators = n_distinct(phenotype_A),
    n_targets   = n_distinct(phenotype_B),
    .groups     = "drop"
  ) %>%
  mutate(significance = "1% FDR")
plot_df <- bind_rows(res_5, res_1) %>%
  left_join(totals, by = "type") %>%
  mutate(
    pct_sig_loc = n_loci / n_loci_tested,
    type        = factor(type, levels = c("forward", "reverse")),
    significance = factor(significance, levels = c("5% FDR", "1% FDR"))
  )

plot_df

### plot

# Split the data by type
plot_df_split <- split(plot_df, plot_df$type) 

# Define your palette (avoid calling it `colors` to not clash with base::colors)
fdr_cols <- c(
  "1% FDR" = "#f39c12",
  "5% FDR" = "#efd4a8" 
)

# Forward plot
df_fwd <- plot_df_split$forward

pr_fwd <- ggplot(df_fwd, aes(x = "", y = pct_sig_loc)) +
  # 5% FDR on top
  geom_col(
    data  = df_fwd %>% dplyr::filter(significance == "5% FDR"),
    aes(fill = significance),
    position = "identity",
    width    = 0.7,
    colour   = "grey40",
    linewidth=0.1
  ) +
  # 1% FDR at the bottom
  geom_col(
    data  = df_fwd %>% dplyr::filter(significance == "1% FDR"),
    aes(fill = significance),
    position = "identity",
    width    = 0.7,
    colour   = "grey40",
    linewidth=0.1
    
  ) +
  geom_text(
    data  = df_fwd %>% dplyr::filter(significance == "5% FDR"),
    aes(label = n_loci),
    vjust    = -0.3,
    size     = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = fdr_cols) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(
    title = "Forward Mediation",
    x     = NULL,
    y     = "% loci with significant triads"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title     = element_text(hjust = 0.5, face = "bold"),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    legend.position = "top"
  )


df_rev <-plot_df_split$rev
pr_rev <- ggplot(df_rev, aes(x = "", y = pct_sig_loc)) +
  # 5% FDR on top
  geom_col(
    data  = df_rev %>% dplyr::filter(significance == "5% FDR"),
    aes(fill = significance),
    position = "identity",
    width    = 0.7,
    colour   = "grey40",
    linewidth=0.1
  ) +
  # 1% FDR at the bottom
  geom_col(
    data  = df_rev %>% dplyr::filter(significance == "1% FDR"),
    aes(fill = significance),
    position = "identity",
    width    = 0.7,
    colour   = "grey40",
    linewidth=0.1
    
  ) +
  geom_text(
    data  = df_rev %>% dplyr::filter(significance == "5% FDR"),
    aes(label = n_loci),
    vjust    = -0.3,
    size     = 5,
    fontface = "bold"
  ) +
  # <-- this is where your colors are applied
  scale_fill_manual(values = fdr_cols) +
  scale_y_continuous(expand = c(0.01, 0)) +
  labs(
    title = "Reverse Mediation",
    x     = NULL,
    y     = "% loci with significant triads"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title     = element_text(hjust = 0.5, face = "bold"),
    axis.text.x    = element_blank(),
    axis.ticks.x   = element_blank(),
    legend.position = "top"
  )



plot_df_fwd <- plot_df_split$forward %>%
  pivot_longer(
    cols = c(n_triads, n_loci, n_mediators, n_targets ),
    names_to = "metric",
    values_to = "value") %>%
  mutate(metric=factor(metric, levels = c("n_triads","n_loci", "n_mediators", "n_targets"))
         
  )
plot_df_rev <- plot_df_split$reverse %>%
  pivot_longer(
    cols = c(n_triads, n_loci, n_mediators, n_targets),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric=factor(metric, levels = c("n_triads","n_loci", "n_mediators", "n_targets"))
         
  )


pr1 <- plot_df_rev %>% filter(significance=="5% FDR") %>%
  filter(metric %in% c( "n_targets", "n_mediators", "n_triads")) %>%
  mutate(metric = factor(metric, levels=c( "n_targets", "n_mediators", "n_triads")))%>%
  ggplot(aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.8, fill="#efd4a8", colour   ="black", linewidth=0.1) +
  geom_text(
    aes(label = scales::comma(value)),      # numeric labels formatted
    vjust = -0.4, size = 5, fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Metric") +
  labs(
    title = "Mediation Analysis (Reverse)",
    y     = NULL,
    x     = NULL
  )+
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_classic() +
  theme(
    axis.title  = element_text(size = 18),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text   = element_text(size = 16),
    axis.text.x =  element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "top"
  )  + coord_flip() 


pr2<-plot_df_fwd %>% filter(significance=="5% FDR") %>%
  filter(metric %in% c( "n_targets", "n_mediators", "n_triads")) %>%
  mutate(metric = factor(metric, levels=c( "n_targets", "n_mediators", "n_triads")))%>%
  ggplot( aes(metric,value, fill =metric)) +
  geom_col(width = 0.8, fill="#efd4a8", colour   ="black", linewidth=0.1) +
  geom_text(
    aes(label = scales::comma(value)),      # numeric labels formatted
    vjust = -0.4, size = 5, fontface = "bold"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Metric") +
  labs(
    title = "Mediation Analysis (Forward)",
    y     = NULL,
    x     = NULL
  ) +
  scale_y_continuous(expand = c(0.01, 0)) +
  theme_classic() +
  theme(
    axis.title  = element_text(size = 18),
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text   = element_text(size = 16),
    axis.text.x =  element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "top"
  ) + coord_flip() 



plt.mediation<-plot_grid(pr_fwd, pr2, pr_rev,pr1,ncol=2, rel_widths = c(0.7, 1), align = "h")
plt.mediation
ggsave("~/cd4_qtl_paper_figures/figure_3/plots/dec2025/mediation_results_summary.pdf", plt.mediation, height = 30, width = 20, units = "cm", dpi = 300)


## 1) Summaries ---------------------------------------------------------

inst_phenotype_A <- findr_res_sig %>%
  dplyr::group_by(type, locus_id) %>%
  dplyr::summarise(
    n_phenotype_A = dplyr::n_distinct(phenotype_A),
    .groups       = "drop"
  )

phenotype_A_to_B <- findr_res_sig %>%
  dplyr::group_by(type, locus_id, phenotype_A) %>%
  dplyr::summarise(
    n_phenotype_B = dplyr::n_distinct(phenotype_B),
    .groups       = "drop"
  )

## 2) Common axis settings ---------------------------------------------

x_limits  <- c(0.9, 10)
x_breaks  <- seq(0, 10, by = 5)


## 3) Helper to make histograms consistent -----------------------------
make_hist <- function(
    df,
    xvar,
    title,
    xlab,
    ylab,
    x_limits  = NULL,
    x_breaks  = NULL,
    fill      = "#716f9e",
    binwidth  = 1,
    line_col  = "#3f3f3f",
    y_limits  = c(0, 110),
    y_breaks  = seq(0, 110, by = 25),
    range_var = NULL   # name of column to show range for (character or NULL)
) {
  # Base plot
  p <- ggplot(df, aes(x = .data[[xvar]])) +
    geom_histogram(
      binwidth  = binwidth,
      fill      = fill,
      colour    = line_col,
      linewidth = 0.2,
      boundary  = 0
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      x     = xlab,
      y     = ylab
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 18),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text  = element_text(size = 16)
    )
  
  # Optional range annotation
  if (!is.null(range_var)) {
    rng <- range(df[[range_var]], na.rm = TRUE)
    label_txt <- paste0("Range: [", round(rng[1], 2), " – ", round(rng[2], 2), "]")
    
    # Use current limits (after scale_y) to place text in top-right
    y_max <- y_limits[2]
    x_max <- if (!is.null(x_limits)) x_limits[2] else max(df[[xvar]], na.rm = TRUE)
    
    p <- p +
      annotate(
        "text",
        x = x_max,
        y = y_max,
        label = label_txt,
        hjust = 1.05,
        vjust = 1.2,
        size  = 4
      )
  }
  
  p
}


## 4) Split data by type -----------------------------------------------

inst_split  <- split(inst_phenotype_A, inst_phenotype_A$type)
pheno_split <- split(phenotype_A_to_B, phenotype_A_to_B$type)

## 5) Build separate graphs (no facets) --------------------------------

# Peaks per Variant (y = # Variants)
p1_forward <- make_hist(
  df       = inst_split$forward,
  xvar     = "n_phenotype_A",
  title    = "Chromatin → Gene",
  xlab     = "No. Unique Peaks",
  ylab     = "No. Unique Variants",
  x_limits = x_limits,
  x_breaks = x_breaks,
  range_var="n_phenotype_A"
)

p1_reverse <- make_hist(
  df       = inst_split$reverse,
  xvar     = "n_phenotype_A",
  title    = "Gene → Chromatin",
  xlab     = "No. Unique Genes",
  ylab     = "No. Unique Variants",
  x_limits = x_limits,
  x_breaks = x_breaks,
  range_var="n_phenotype_A"
  
)

# Genes per Peak (y = # Peaks)
p2_forward <- make_hist(
  df       = pheno_split$forward,
  xvar     = "n_phenotype_B",
  title    = "Gene → Chromatin",
  xlab     = "No. Unique Genes",
  ylab     = "No. Unique Peaks",
  x_limits = x_limits,
  x_breaks = x_breaks,
  range_var="n_phenotype_B"
  
)

p2_reverse <- make_hist(
  df       = pheno_split$reverse,
  xvar     = "n_phenotype_B",
  title    = "Chromatin → Gene",
  xlab     = "No. Unique Peaks",
  ylab     = "No. Unique Genes",
  x_limits = x_limits,
  x_breaks = x_breaks,
  range_var="n_phenotype_B"
  
)

## 6) Grid layout (2 × 2) ----------------------------------------------

combined_plot <- cowplot::plot_grid(
  p1_forward, p2_forward,
  p1_reverse, p2_reverse,
  ncol  = 2,
  align = "hv"
)

combined_plot


## 4. Grid layout (side‑by‑side) --------------------------------
combined_plot <- plot_grid(
  p1_forward, p2_forward, p1_reverse, p2_reverse,
  ncol        = 2,   # one column → vertical stack
  align       = "v", # align plots vertically
  axis        = "lr",  # keep left/right axes lined up
  rel_heights = c(1,1,1,1)
)
combined_plot


combined_plot<-plot_grid(pr_fwd, pr2, p1_forward, p2_forward,
                         pr_rev,pr1, p1_reverse, p2_reverse,
                         nrow=2, rel_widths = c(.6, 1, 0.8, 0.8), align = "h")

ggsave("~/cd4_qtl_paper_figures/figure_3/plots/dec2025/causal_inference_summary_distributions_grid.pdf", combined_plot, width = 32, height = 20, units = "cm", dpi = 300)


library(VennDiagram)
library(grid)
# Per-instrument/type counts (you already have this as `n`)
counts <- findr_res_sig %>% mutate(triad=ifelse(type == "forward", paste0(locus_id, "_", phenotype_A, "_", phenotype_B), paste0(locus_id, "_", phenotype_B, "_", phenotype_A))) %>%
  dplyr::count(triad, type, name = "n")

forward_ids <- counts %>% filter(type == "forward") %>% pull(triad) %>% unique()

reverse_ids <- counts %>% filter(type == "reverse") %>% pull(triad) %>% unique()

grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1      = length(forward_ids),
  area2      = length(reverse_ids),
  cross.area = length(intersect(forward_ids, reverse_ids)),
  category   = c("Forward", "Reverse"),
  fill       = c("#56B4E9", "#E69F00"),
  alpha      = 0.5,
  cex        = 1.3,
  cat.cex    = 1.2
)

pdf("~/cd4_qtl_paper_figures/figure_3/plots/dec2025/mediation_FDR5_results_venn_instrument_direction.pdf")
venn<-grid.draw(venn.plot)
dev.off()

