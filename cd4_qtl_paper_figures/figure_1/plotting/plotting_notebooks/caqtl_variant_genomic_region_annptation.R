.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4", "/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(GenomicRanges)
  library(VariantAnnotation)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(cowplot)
  library(readr)
  library(stringr)
  library(RColorBrewer)
})

source("~/cd4_qtl_paper_figures/utils/color_pallete_helper.R")

summary_df<-fread("/gchm/cd4_qtl_paper_figures/figure_1/data/eqtl_caqtl_finemapping_coloc_all.tsv")

## one filter, reused everywhere
caqtl_only <- c("coloc", "caQTL_only")
caqtl_df0 <- summary_df %>% mutate(caqtl_category=ifelse(cs_type_any_var %in% c("in_other_Peak","in_corr_Peak"), "in_other_Peak", cs_type_any_var)) %>%
  filter(coloc_status %in% caqtl_only) 
rm(summary_df)
variantsGR.caqtl <- GRanges(
  seqnames = caqtl_df0$chr.x,
  ranges   = IRanges(start = caqtl_df0$variant_pos.x, end = caqtl_df0$variant_pos.x)
)

mcols(variantsGR.caqtl)$finemapped_cs_caqtl <- caqtl_df0$finemapped_cs_caqtl
mcols(variantsGR.caqtl)$caqtl_category     <- caqtl_df0$caqtl_category
mcols(variantsGR.caqtl)$variant_id.x    <- caqtl_df0$variant_id.x


# ======================================================================
# 3) GENOMIC REGION (TxDb locateVariants → CS-level majority vote)
# ======================================================================
# PRO TIP: build the TxDb once, cache to SQLite, and reload in future runs.
txdb <- GenomicFeatures::makeTxDbFromGFF(
  "~/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf"
)

## --- 0) hygiene: make seqlevels consistent & keep standard chromosomes
dups <- duplicated(variantsGR.caqtl)  # compares seqnames/start/end/strand
variantsGR.caqtl_unique <- variantsGR.caqtl[!dups]

seqlevelsStyle(variantsGR.caqtl_unique) <- seqlevelsStyle(txdb)
variantsGR.caqtl_unique <- keepStandardChromosomes(variantsGR.caqtl_unique, pruning.mode = "coarse")
variantsGR.caqtl_unique <- sort(variantsGR.caqtl_unique)

## --- 1) define the variant classes you want to locate
class_map <- list(
  coding     = CodingVariants(),
  fiveUTR    = FiveUTRVariants(),
  threeUTR   = ThreeUTRVariants(),
  spliceSite = SpliceSiteVariants(),
  intron     = IntronVariants(),
  promoter   = PromoterVariants(upstream = 2000, downstream = 200),
  intergenic = IntergenicVariants()
)

## --- 2) chunk the query GRanges
chunk_size <- 10000L  # tune to your memory; smaller = safer
cuts       <- ceiling(seq_along(variantsGR.caqtl_unique) / chunk_size)
gr_chunks  <- split(variantsGR.caqtl_unique, cuts)

## --- 3) helper that runs locateVariants for all classes on one chunk
run_locate_all <- function(qgr) {
  lapply(class_map, function(cls) {
    locateVariants(qgr, txdb, cls)
  })
}

## --- 4A) SERIAL loop (most robust)
hits_list <- lapply(names(class_map), function(nm) GRanges())
names(hits_list) <- names(class_map)

for (i in seq_along(gr_chunks)) {
  qgr <- gr_chunks[[i]]
  message(sprintf("Chunk %d/%d (n=%d variants)", i, length(gr_chunks), length(qgr)))
  chunk_res <- tryCatch(run_locate_all(qgr),
                        error = function(e) { message("  -> ERROR: ", conditionMessage(e)); NULL })
  if (is.null(chunk_res)) next
  # append results per class
  for (nm in names(class_map)) {
    hits_list[[nm]] <- c(hits_list[[nm]], chunk_res[[nm]])
  }
}



region_priority <- c("coding","fiveUTR","threeUTR","spliceSite","promoter","intron","intergenic")

hits_tbl <- rbindlist(lapply(names(hits_list), function(lbl) {
  hmc <- as.data.frame(S4Vectors::mcols(hits_list[[lbl]]))
  if (!"QUERYID" %in% names(hmc)) stop("locateVariants result missing QUERYID")
  data.frame(
    QUERYID = as.integer(hmc$QUERYID),
    TXID    = if ("TXID" %in% names(hmc)) as.character(hmc$TXID) else NA_character_,
    GENEID  = if ("GENEID" %in% names(hmc)) as.character(hmc$GENEID) else NA_character_,
    region  = lbl,
    stringsAsFactors = FALSE
  )
}), use.names = TRUE)

saveRDS(hits_tbl, "~/cd4_qtl_paper_figures/figure_1/data/caqtl_variant_genomic_annotation.rds")