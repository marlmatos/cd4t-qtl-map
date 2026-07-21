# ============================================================
# author: Marliette Matos
# Annotate chromatin-accessibility peaks by genomic context and
# ChromBPNet motif content.
#
# The script:
#   1. Loads the filtered set of chromatin-accessibility peaks.
#   2. Classifies each peak as promoter or distal based on
#      direct overlap with a gene transcription start site.
#   3. Imports ChromBPNet motif instances from a BED file.
#   4. Counts the number of occurrences of each motif within
#      every accessibility peak.
#   5. Constructs a peak-by-motif count matrix and appends the
#      motif counts to the peak annotation table.
#   6. Saves the resulting peak-level annotations for use in
#      downstream enrichment and visualization analyses.
#
# Promoter status is defined using exact TSS overlap without an
# additional flanking window; all remaining peaks are classified
# as distal regulatory elements.
# ============================================================


source("~/cd4t-qtl-map/analysis/paths_config.R")

.libPaths(fig1_paths$r_env)

library(ggplot2)
library(dplyr)
library(readr) 
library(data.table)
library(GenomicRanges)
library(tibble)
library(forcats)

options(bitmapType = "cairo")

peak_coords <- read_tsv(fig1_paths$peak_coords_file, col_names = NULL)
peaks_names<-read.delim2(fig1_paths$peaks_names_file)$peak_name
peak_coords <-peak_coords[peak_coords$X4 %in% peaks_names,]
peak_coordsGR <- GRanges(peak_coords$X1, IRanges(peak_coords$X2, peak_coords$X3), peak_name = peak_coords$X4)


names(peak_coordsGR) <- peak_coordsGR$peak_name

#read the gencode gene annotation
gencode <- fread(fig1_paths$gencode_file,
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
ov_tss  <- countOverlaps(peak_coordsGR, tss_gr, ignore.strand = TRUE)


peak_annotation <- tibble(
  peak_name = names(peak_coordsGR),
  chr       = as.character(seqnames(peak_coordsGR)),
  start     = start(peak_coordsGR),
  end       = end(peak_coordsGR),
  promoter  = ov_tss > 0,
  annotation = if_else(promoter, "promoter", "distal")
)


# peek
dplyr::count(peak_annotation, annotation)
saveRDS(peak_annotation, fig1_paths$peak_annotation_file)

saveRDS(peak_annotation, fig1_paths$peak_annotation_no_window_file)


#### peak overlaps with 
library(rtracklayer)
library(GenomicRanges)
library(Matrix)
library(dplyr)

chrombpnet_motifs <- rtracklayer::import(fig1_paths$chrombpnet_motifs_file, format="BED")

ov_motif  <- countOverlaps(peak_coordsGR, chrombpnet_motifs, ignore.strand = TRUE)

# indices
i <- queryHits(ov_motif)        # rows: peaks
s <- subjectHits(ov_motif)      # index into chrombpnet_motifs

# column indices by motif name (factor encodes columns)
motif_names <- mcols(chrombpnet_motifs)$name
motif_fac   <- factor(motif_names)              # unique columns per motif
j <- as.integer(motif_fac[s])                   # columns: motifs (via subject indices)

# dims
n_peaks  <- length(peak_coordsGR)
n_motifs <- nlevels(motif_fac)

motif_mat <- sparseMatrix(
  i    = i,
  j    = j,
  x    = 1L,
  dims = c(n_peaks, n_motifs),
  dimnames = list(
    NULL,
    levels(motif_fac)
  )
)

motif_mat

peak_name <- names(peak_coordsGR)
if (is.null(peak_name)) {
  peak_name <- paste0(seqnames(peak_coordsGR), ":", start(peak_coordsGR), "-", end(peak_coordsGR))
}

stopifnot(nrow(motif_mat) == length(peak_name))

## 4A)Dense wide table then join
motif_counts <- as_tibble(as.matrix(motif_mat)) |>
  add_column(peak_name = peak_name, .before = 1)

peak_annotation <- peak_annotation |>
  left_join(motif_counts, by = "peak_name")
# Replace NAs in motif columns with 0 if desired:
motif_cols <- setdiff(names(motif_counts), "peak_name")
peak_annotation <- peak_annotation |>
  mutate(across(all_of(motif_cols), ~ tidyr::replace_na(., 0L)))
saveRDS(peak_annotation, fig1_paths$peak_annotation_no_window_file)
