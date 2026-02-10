library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)
library(GenomicRanges)
library(cowplot)

options(bitmapType = "cairo")
source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")
#read all finemapped variants
summary_df<-read_delim("/gpfs/commons/home/mmatos/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")
#read chromBPNet results annoated with motifs hits
cbpnet_variants<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/cd4_top_cpbnet_variants_motif_overlap_moreTFS_boolean.tsv")
#read variant_motif_intercept
variant_motif_map<-read_tsv("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/motif_variant_overlap_all.tsv")


### what fraction of chromBPNet variants are eQTL/caQTls/both?
# get the chromBPNet variants that are within caQTLs credible sets
caqtl_fract<-summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.x)
  ) %>%
  semi_join(cbpnet_variants, by = c("var_id" = "variant_id")) %>%
  group_by(coloc_status) %>%
  summarise(var_ids = list(unique(var_id)), .groups = "drop")

# get the chromBPNet variants that are within eQTLs credible sets
eqtl_fract<-summary_df %>%
  mutate(
    var_id = gsub("([0-9]+):([0-9]+)\\[b38\\]([A-Z]),([A-Z])",
                  "chr\\1_\\2_\\3_\\4",
                  variant_id.y)
  ) %>%
  semi_join(cbpnet_variants, by = c("var_id" = "variant_id")) %>%
  group_by(coloc_status)  %>%
  summarise(var_ids = list(unique(var_id)), .groups = "drop")

rm(summary_df)

##----------------
#intercept both lists to get unique "coloc" variants
# get shared coloc variants
ca_coloc <- caqtl_fract %>%
  filter(coloc_status == "coloc") %>%
  pull(var_ids) %>% .[[1]]

eq_coloc <- eqtl_fract %>%
  filter(coloc_status == "coloc") %>%
  pull(var_ids) %>% .[[1]]

shared_coloc <- intersect(ca_coloc, eq_coloc)
length(shared_coloc)         # how many


# get unique fractions
ca_only <- caqtl_fract %>%
  dplyr::filter(coloc_status == "caQTL_only") %>%
  dplyr::pull(var_ids) %>% .[[1]]

eq_only <- eqtl_fract %>%
  dplyr::filter(coloc_status == "eQTL_only") %>%
  dplyr::pull(var_ids) %>% .[[1]]

##make a composite list
res <- list(
  ca_only   = unique(ca_only),
  eq_only   = unique(eq_only),
  shared    = unique(shared_coloc)
)

# quick sanity counts
sapply(res, length)

#sanity check to make sure these are muatuallu exclusive variant lists
# 1) Pairwise overlaps must be zero
all(combn(names(res), 2, \(p) length(intersect(res[[p[1]]], res[[p[2]]]))) == 0)

# 2) Union size should equal sum of sizes
length(unique(unlist(res, use.names = FALSE))) == sum(sapply(res, length))


#### depuplicate if shared
#move anything that appears in more than one list into shared
#remove those from ca_only and eq_only
#deduplicate each bucket


dedup_buckets <- function(res) {
  ca0 <- unique(na.omit(res$ca_only))
  eq0 <- unique(na.omit(res$eq_only))
  sh0 <- unique(na.omit(res$shared))
  
  stacked <- tibble(
    var_id = c(ca0, eq0, sh0),
    source = c(rep("caQTL_only", length(ca0)),
               rep("eQTL_only", length(eq0)),
               rep("shared", length(sh0)))
  )
  
  # vars present in >= 2 lists
  multi <- stacked %>%
    distinct(var_id, source) %>%
    add_count(var_id, name = "lists_containing") %>%
    filter(lists_containing >= 2) %>%
    pull(var_id) %>%
    unique()
  
  shared_new  <- unique(c(sh0, multi))
  ca_only_new <- setdiff(ca0, shared_new)
  eq_only_new <- setdiff(eq0, shared_new)
  
  list(
    caQTL_only = ca_only_new,
    eQTL_only = eq_only_new,
    shared  = shared_new
  )
}

res_dedup <- dedup_buckets(res)

# sanity checks
all(combn(names(res_dedup), 2, \(p) length(intersect(res_dedup[[p[1]]], res_dedup[[p[2]]]))) == 0)
length(unique(unlist(res_dedup, use.names = FALSE))) == sum(sapply(res_dedup, length))
sapply(res_dedup, length)

##add the chrombpnet only variants to the list
already_bucketed <- unique(c(res_dedup$caQTL_only,
                             res_dedup$eQTL_only,
                             res_dedup$shared))
cbpnet_all <- unique(na.omit(cbpnet_variants$variant_id))

# 3) Variants seen in CBPNet but in none of the buckets
cbpnet_only <- setdiff(cbpnet_all, already_bucketed)

# 4) Add to the list
res_dedup$chromBPnet <- cbpnet_only

##get summary counts for overlap with QTLs finemapping
sapply(res_dedup, length)

#join all results in tibble df
cbpnet_vars_qtl_annot <- bind_rows(
  lapply(res_dedup, \(v) tibble(var_id = v)),
  .id = "bucket"
)

#join the motif overls
cbpnet_vars_qtl_annot <- cbpnet_vars_qtl_annot %>% left_join(variant_motif_map, join_by(var_id==variant_id)) #notice that the df lenght went up by 56

table(duplicated(cbpnet_vars_qtl_annot$var_id)) #56 duplicated vars

#56 variants have more than one overlapping motif
mutiple_hits<-cbpnet_vars_qtl_annot %>%
  group_by(var_id) %>%
  summarise(
    n_overlaps       = n(),                        # total rows (overlaps)
    n_unique_motifs  = n_distinct(TF_family),      # or use motif_id if you have it
    has_multiple     = n_unique_motifs > 1,
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_motifs))


##difference in QTL breakdown if motif/nomotif
cbpnet_vars_qtl_annot<-cbpnet_vars_qtl_annot %>%
  mutate(motif_disrupt = !is.na(motif_id),
         motif_disrupt = factor(motif_disrupt, levels = c(FALSE, TRUE),
                                labels = c("No Motif", "Motif Overlap")))



plot_df <- cbpnet_vars_qtl_annot %>%
  group_by(motif_disrupt, bucket) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(bucket = factor(bucket, levels = c("chromBPnet","shared","caQTL_only","eQTL_only")))

# counts, side-by-side
pt1<-ggplot(plot_df, aes(x = motif_disrupt, y = count, fill = bucket)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +  # <- stat="identity"
  labs(x = NULL, y = "Variants", fill = "Variant set") +
  theme_classic() + theme(legend.position = "none")

# proportions within each x (optional)
pt2<-ggplot(plot_df, aes(x = motif_disrupt, y = count, fill = bucket)) +
  geom_col(position = "fill", width = 0.75) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of variants", fill = "Variant set") +
  theme_classic()

grid<-plot_grid(pt1,pt2, rel_widths = c(1,0.5), nrow = 1,
          align = "h",
          axis = "tb")
ggsave("~/cd4_qtl_paper_figures/figure_2_chrombpnet/plotting/variant_motif_overlap/plots/variant_motif_overlap_qtl.pdf", grid, width=10, height = 4)

print("Variant-TF motif Overlap:")
table(cbpnet_vars_qtl_annot$motif_disrupt)


#### Plot breakdown per motif family

sum_df<-cbpnet_vars_qtl_annot %>%
  filter(motif_disrupt=="Motif Overlap") %>%
  group_by(TF_family, bucket) %>%
  summarise(n_vars = n_distinct(var_id), .groups = "drop") %>%
  group_by(TF_family) %>%
  mutate(total = sum(n_vars)) %>%            # for ordering families by overall size
  ungroup() %>%
  mutate(TF_family = fct_reorder(TF_family, n_vars),
         bucket = factor(bucket,
                         levels = c( "chromBPnet", "shared", "caQTL_only", "eQTL_only" )))  # sort by total


motif_qtl.plt<-ggplot(sum_df, aes(x = TF_family, y = n_vars, fill = bucket)) +
  geom_col(position = "fill", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = variant_cmap2, drop = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of motif overlapping variants", fill = "Variant Set") +
  theme_classic()


library(ggpubr)

cbpnet_vars_qtl_annot %>%
  filter(motif_disrupt == "Motif Overlap") %>%
  mutate(bucket2 = ifelse(bucket != "chromBPnet", "QTL", bucket)) %>%
  ggplot(aes(y = hit_coefficient_global, x = bucket2)) +
  geom_boxplot() +
  stat_compare_means(
    method = "t.test",
    label = "p.format",   # formats p-value nicely
    comparisons = list(c("chromBPnet", "QTL"))
  ) +
  labs(x = NULL, y = "Hit coefficient (global)")



#########--- check for motif overlap on peak annotations (rpomoter vs. distal)

#read atac called peaks for chromBPnet  
bp.peaks<-fread("cd4_chrombpnet/data/inputs/peaks/merged_cd4_samples_peaks_no_blacklist.sorted.bed")
bp.peaksGR <- GRanges(bp.peaks$V1, IRanges(bp.peaks$V2, bp.peaks$V3), peak_name = bp.peaks$V4)

#read the gencode gene annotation
gencode <- fread("~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.genes.bed",
                 sep="\t", data.table=FALSE) %>%
  #filter(V8 == "protein_coding",
  #      grepl("^chr(\\d+|X|Y)$", V1)) %>%                 # drop alts/patches
  transmute(
    chr        = V1,
    strand     = V6,
    gene_id    = V7,
    gene       = V4,                                       # gene symbol
    tss_1based = if_else(V6 == "+", V2 + 1L, V3)           # BED start is 0-based
  ) %>%
  distinct(chr, gene, .keep_all = TRUE)                    # ensure 1 row per chr×gene


###### data driven annotation of regulatory elements ####
### promoter: peak needs to overlap at least one TSS of a gene
### distal regulatory element -- no TSS overlap


# 1) Build 1-bp TSS GRanges
tss_gr <- with(gencode,
               GRanges(
                 seqnames = chr,
                 ranges   = IRanges(start = tss_1based, end = tss_1based),
                 strand   = strand
               )
)

# 5) Overlap: promoter if overlaps ANY TSS (ignore strand for promoter call)
#tss_win <- resize(tss_gr, width = 201, fix = "center")  # ±100bp
ov_tss  <- countOverlaps(bp.peaksGR, tss_gr, ignore.strand = TRUE)

# Write annotations back onto the GRanges
mcols(bp.peaksGR)$promoter   <- ov_tss > 0
mcols(bp.peaksGR)$annotation <- ifelse(mcols(bp.peaksGR)$promoter, "promoter", "distal")

## now annotate motif_vars
cbpnet_variants_GR<-readRDS("~/cd4_qtl_paper_figures/figure_2_chrombpnet/data/motif_variant_overlap/cd4_top_cpbnet_variants_granges_motif_overlap_moreTFS_boolean.rds")

# 1) Overlap mapping (variant -> peak)
hits <- findOverlaps(cbpnet_variants_GR, bp.peaksGR, ignore.strand = TRUE)

## 3) Transfer peak metadata to variants
# Prepare vectors aligned to variants
ann_vec  <- rep(NA_character_, length(cbpnet_variants_GR))
name_vec <- rep(NA_character_, length(cbpnet_variants_GR))

ann_vec[queryHits(hits)]  <- mcols(bp.peaksGR)$annotation[subjectHits(hits)]
name_vec[queryHits(hits)] <- mcols(bp.peaksGR)$peak_name[subjectHits(hits)]

mcols(cbpnet_variants_GR)$peak_annotation <- ann_vec
mcols(cbpnet_variants_GR)$peak_name       <- name_vec


var_tss <- tibble(
  var_id= cbpnet_variants_GR$variant_id,
  peak_annotation=cbpnet_variants_GR$peak_annotation
)

##join with back to the chrombpnet var df contaitning the motifs ids
cbpnet_vars_qtl_annot <- cbpnet_vars_qtl_annot %>%  left_join(variant_annot_df, by = "var_id")


sum_df2 <- cbpnet_vars_qtl_annot %>%
  filter(!is.na(TF_family), !is.na(peak_annotation)) %>%
  group_by(TF_family,peak_annotation) %>%
  summarise(n_vars = n_distinct(var_id), .groups = "drop_last") %>%
  # order families by their total across buckets/annotations
  group_by(TF_family) %>%
  mutate(total = sum(n_vars)) %>%
  ungroup() %>%
  mutate(
    TF_family = fct_reorder(TF_family, total))
  

# 3) Plot: 100% stacked by peak annotation, faceted by bucket
motif_tss.plt <- ggplot(sum_df2, aes(x = TF_family, y = n_vars, fill = peak_annotation)) +
  geom_col(position = "fill", width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = "Proportion of motif-overlapping variants",
    fill = "Genomic annotation"
  ) +
  theme_classic()

motif_tss.plt

