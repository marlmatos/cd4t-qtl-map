#!/usr/bin/env Rscript
.libPaths(c("/gchm/R/x86_64-pc-linux-gnu-library/4.4","/nfs/sw/easybuild/software/R/4.4.1-gfbf-2023b/lib/R/library"))

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})


in_path  <- "/gchm/cd4_chrombpnet/chrombpnet_model_b7/variant_prediction_scores/averaged_scores/average_cd4_tcells_AJ_common_variants.mean.variant_scores_IPS05.tsv"

# dbSNP beds (same style/locations as your SuSiE script)
dbsnp38_dir <- "/gcgl/sghatan/genome/dbsnp151_GRCh38p7"   # bed_chr_1.bed.gz ... bed_chr_22.bed.gz (+ _X)
dbsnp37_dir <- "/gcgl/sghatan/genome/dbsnp151_GRCh37p13"  # bed_chr_1.bed.gz ... bed_chr_22.bed.gz (+ _X)

# Where to drop intermediates + outputs  (use absolute path; expand ~)
base_out <- path.expand("~/cd4_qtl_paper_figures/figure_6/data/liftover_cbpnet_rsID")
bed_dir  <- file.path(base_out, "bed38")
int_dir  <- file.path(base_out, "intersect38")
out_dir  <- file.path(base_out, "final")
dir.create(bed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Path to bedtools (if not on PATH, edit)
bedtools <- Sys.getenv("BEDTOOLS", unset="/gchm/.conda/envs/liftover/bin/bedtools")

# Sanity: ensure bedtools exists
if (!file.exists(bedtools)) stop("bedtools not found at: ", bedtools)

# ------------------ 1) Read cbpnet (your way) ------------------
cbpnet <- fread(in_path, data.table = FALSE, stringsAsFactors = FALSE) %>%
  dplyr::select(variant_id)
stopifnot("variant_id" %in% names(cbpnet))

# ------------------ 2) Parse variant_id -> chr/pos/ref/alt ------------------
parts <- str_split_fixed(cbpnet$variant_id, "_", 4)
if (ncol(parts) < 4) stop("variant_id not in 'chr_pos_ref_alt' format")

cbpnet$chr <- parts[,1]
cbpnet$pos <- as.integer(parts[,2])
cbpnet$ref <- parts[,3]
cbpnet$alt <- parts[,4]

# Ensure chr has "chr" prefix
cbpnet$chr <- ifelse(startsWith(cbpnet$chr, "chr"), cbpnet$chr, paste0("chr", cbpnet$chr))

# Limit to autosomes + X present in file
chr_levels <- intersect(unique(cbpnet$chr), c(paste0("chr", 1:22), "chrX"))
if (length(chr_levels) == 0L) stop("No autosomes/X found in variant_id column.")

# ------------------ 3) Build per-chr BEDs (hg38) ------------------
# BED columns: chr, start(0-based), end(1-based), id, ref, alt
write_bed <- function(df_chr, out_bed) {
  bed <- data.frame(
    chr   = df_chr$chr,
    start = as.integer(df_chr$pos - 1L),
    end   = as.integer(df_chr$pos),
    id    = df_chr$variant_id,
    ref   = df_chr$ref,
    alt   = df_chr$alt,
    stringsAsFactors = FALSE
  )
  fwrite(bed, out_bed, sep = "\t", quote = FALSE, col.names = FALSE)
}

by_chr <- split(cbpnet, cbpnet$chr)
for (c in chr_levels) {
  this_bed <- file.path(bed_dir, paste0(c, ".hg38.bed"))
  if (nrow(by_chr[[c]]) > 0) write_bed(by_chr[[c]], this_bed)
}

# ------------------ 4) Intersect with dbSNP GRCh38p7 to get rsIDs ------------------
# Use system2() to avoid shell redirection/tilde issues.
intersections <- list()

for (c in chr_levels) {
  bed_a  <- file.path(bed_dir, paste0(c, ".hg38.bed"))
  if (!file.exists(bed_a) || file.size(bed_a) == 0) next
  
  chr_tag <- sub("^chr", "", c)  # "chr10" -> "10"
  bed_b   <- file.path(dbsnp38_dir, paste0("sorted_bed_chr_", chr_tag, ".bed.gz"))
  out_txt <- file.path(int_dir, paste0(c, ".rs38.txt"))
  
  status <- system2(
    bedtools,
    args   = c("intersect", "-wa", "-wb", "-a", bed_a, "-b", bed_b),
    stdout = out_txt
  )
  
  if (status != 0L) {
    warning("bedtools intersect failed for ", c, " with status ", status)
  }
  
  if (file.exists(out_txt) && file.size(out_txt) > 0) {
    dt <- fread(out_txt, header = FALSE)
    # our 6 cols + dbSNP 6 cols
    setnames(dt, c("a_chr","a_start","a_end","a_id","a_ref","a_alt",
                   "b_chr","b_start","b_end","rsid","score","strand"))
    intersections[[c]] <- dt
  }
}


rs38_all <- if (length(intersections)) rbindlist(intersections, use.names = TRUE, fill = TRUE) else data.table()

# ------------------ 5) Join rsIDs to dbSNP GRCh37 to get hg19 positions ------------------
if (nrow(rs38_all) == 0L) {
  warning("No variants overlapped dbSNP GRCh38p7; no rsIDs to carry to hg19.")
  res <- cbpnet %>%
    transmute(
      variant_id,
      hg38_chr = chr,
      hg38_pos = pos,
      ref, alt,
      rsid = NA_character_,
      hg19_chr = NA_character_,
      hg19_pos = NA_integer_    )
  fwrite(res, file.path(out_dir, "cbpnet_hg38_to_hg19_rsIDjoin.tsv.gz"), sep = "\t")
  quit(save = "no")
}

# one row per variant_id (if multiple rsIDs overlapped, keep the first by genomic order)
rs38_unique <- rs38_all %>%
  arrange(a_id, b_chr, b_start) %>%
  group_by(a_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(variant_id = a_id,
         hg38_chr = a_chr,
         hg38_pos = a_end,  # 1-based position
         ref = a_ref,
         alt = a_alt,
         rsid)

# Load dbSNP37 per chr we saw; skip track header line like your earlier script
load_dbsnp37_chr <- function(chr_tag) {
  f <- file.path(dbsnp37_dir, paste0("bed_chr_", chr_tag, ".bed.gz"))
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f, skip = 1), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  setnames(dt, c("chr","start.grc37","end.grc37","rsid","score","strand"))
  dt[, .(chr, start.grc37, end.grc37, rsid)]
}

dbsnp37_list <- lapply(sub("^chr", "", chr_levels), load_dbsnp37_chr)
dbsnp37_all  <- rbindlist(Filter(Negate(is.null), dbsnp37_list), use.names = TRUE, fill = TRUE)

# Join by rsid to get hg19 position
res <- rs38_unique %>%
  left_join(dbsnp37_all, by = "rsid") %>%
  transmute(
    variant_id,
    hg38_chr,
    hg38_pos = as.integer(hg38_pos),
    ref, alt,
    rsid,
    hg19_chr = chr,
    hg19_pos = as.integer(end.grc37)
  )

# Merge back cbpnet metrics
res <- cbpnet %>%
  select(variant_id) %>%
  right_join(res, by = "variant_id") %>%
  relocate(variant_id, hg38_chr, hg38_pos, ref, alt, rsid, hg19_chr, hg19_pos)

# ------------------ 6) Also report variants with no rsID / no hg19 ------------------
has_rs   <- !is.na(res$rsid) & res$rsid != ""
has_hg19 <- has_rs & !is.na(res$hg19_chr) & !is.na(res$hg19_pos)

unmapped <- cbpnet %>%
  filter(!(variant_id %in% res$variant_id[has_rs])) %>%
  transmute(variant_id, chr, pos, ref, alt)
job
no_hg19  <- res %>%
  filter(!has_hg19) %>%
  select(variant_id, hg38_chr, hg38_pos, ref, alt, rsid)

# ------------------ 7) Write outputs ------------------
fwrite(res,      file.path(out_dir, "cbpnet_hg38_to_hg19_rsIDjoin.tsv.gz"), sep = "\t")
fwrite(unmapped, file.path(out_dir, "cbpnet_no_dbSNP_rsid.tsv.gz"),        sep = "\t")
fwrite(no_hg19,  file.path(out_dir, "cbpnet_rsid_no_hg19.tsv.gz"),         sep = "\t")

message("✓ Wrote: ", file.path(out_dir, "cbpnet_hg38_to_hg19_rsIDjoin.tsv.gz"))
message("  - No dbSNP rsID for: ", nrow(unmapped))
message("  - rsID present but no hg19 match: ", nrow(no_hg19))
