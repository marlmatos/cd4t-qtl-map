# --- deps ---
library(rhdf5)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(GenomeInfoDb)

# ---------- H5 helpers ----------
.h5_normkey <- function(dset) if (!startsWith(dset, "/")) paste0("/", dset) else dset

.h5_find <- function(h5_path, dset) {
  key <- .h5_normkey(dset)
  tbl <- h5ls(h5_path, recursive = TRUE)
  hit <- tbl[file.path(tbl$group, tbl$name) == key, , drop = FALSE]
  if (nrow(hit) == 0) stop("Dataset not found: ", key)
  hit
}

.h5_dims <- function(h5_path, dset) {
  row <- .h5_find(h5_path, dset)
  if (!nzchar(row$dim[1])) stop("Dataset has no dims? ", dset)
  as.integer(strsplit(row$dim[1], " x ")[[1]])
}

# Read a 2D "profiles" dataset for a single variant; handle (L,N) or (N,L) gracefully.
# Returns a numeric vector of length L (always the genomic-length axis).
.h5_read_profile_1variant <- function(h5_path, dset, idx) {
  key  <- .h5_normkey(dset)
  dims <- .h5_dims(h5_path, key)
  if (length(dims) != 2L) stop("Expected 2D dataset for ", key, "; got dims ", paste(dims, collapse = "x"))
  # First try column selection (L,N) -> [:, idx]
  out <- try({
    L <- dims[1]; N <- dims[2]
    if (idx < 1L || idx > N) stop(sprintf("idx=%d out of range [1,%d] for %s", idx, N, key))
    as.numeric(drop(h5read(h5_path, key, index = list(1:L, as.integer(idx)))))
  }, silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # Fallback: row selection (N,L) -> [idx, :]
  out <- try({
    N <- dims[1]; L <- dims[2]
    if (idx < 1L || idx > N) stop(sprintf("idx=%d out of range [1,%d] for %s (row-major)", idx, N, key))
    as.numeric(drop(h5read(h5_path, key, index = list(as.integer(idx), 1:L))))
  }, silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  stop("Could not read a 1-variant profile from ", key, " with idx=", idx,
       ". Neither (L,N) column nor (N,L) row access worked.")
}

# Flexible scalar read for counts datasets with shapes (1,N), (N,1), (N), or a 2D row/col.
.h5_read_count_scalar <- function(h5_path, dset, idx) {
  key <- .h5_normkey(dset)
  dims <- .h5_dims(h5_path, key)
  if (length(dims) == 1L) {
    N <- dims[1]
    if (idx < 1L || idx > N) stop(sprintf("idx=%d out of range [1,%d] for %s", idx, N, key))
    return(as.numeric(h5read(h5_path, key, index = list(as.integer(idx)))))
  } else if (length(dims) == 2L) {
    # Try (1,N) or (N,1)
    out <- try(as.numeric(h5read(h5_path, key, index = list(1L, as.integer(idx)))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    out <- try(as.numeric(h5read(h5_path, key, index = list(as.integer(idx), 1L))), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    # As a last resort, accept per-position counts and sum them
    # (some pipelines store total counts as sum over profile)
    # Try column:
    out <- try({
      L <- dims[1]; N <- dims[2]
      if (idx >= 1L && idx <= N) {
        vec <- as.numeric(drop(h5read(h5_path, key, index = list(1:L, as.integer(idx)))))
        sum(vec)
      } else stop("bad idx")
    }, silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    # Try row:
    out <- try({
      N <- dims[1]; L <- dims[2]
      if (idx >= 1L && idx <= N) {
        vec <- as.numeric(drop(h5read(h5_path, key, index = list(as.integer(idx), 1:L))))
        sum(vec)
      } else stop("bad idx")
    }, silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
  }
  stop(sprintf("Unable to read count scalar from '%s' at idx=%s; dataset shape=%s not supported.",
               key, idx, paste(dims, collapse="x")))
}

# Load variant_ids from an H5 dataset or text file (1 per line or TSV first col)
.load_variant_ids <- function(src_path, ids_key = "variant_ids") {
  if (grepl("\\.h5$|\\.hdf5$", src_path, ignore.case = TRUE)) {
    key <- ids_key
    if (!H5Lexists(src_path, key)) {
      alt <- if (startsWith(key, "/")) substring(key, 2L) else paste0("/", key)
      if (!H5Lexists(src_path, alt)) stop(sprintf("Dataset '%s' not found in %s", ids_key, src_path))
      key <- alt
    }
    ids <- h5read(src_path, key)
    return(as.character(ids))
  } else {
    con <- file(src_path, "r"); on.exit(close(con), add = TRUE)
    lines <- readLines(con)
    ids <- sub("\t.*$", "", lines)
    ids[nzchar(ids)]
  }
}

# ---------- Public: parse one variant from averaged predictions H5 ----------
# You can address by idx (1-based) or by variant_id (using variant_ids or variant_ids_path).
parse_variant_from_h5 <- function(
    h5_file_path,
    idx = NULL,
    variant_id = NULL,
    variant_ids = NULL,
    variant_ids_path = NULL,
    ids_key = "variant_ids",
    group = "observed"
) {
  stopifnot(file.exists(h5_file_path))
  if (is.null(idx) && is.null(variant_id)) stop("Provide either idx or variant_id.")
  if (!is.null(variant_id) && is.null(idx)) {
    if (is.null(variant_ids)) {
      if (is.null(variant_ids_path)) stop("Provide variant_ids or variant_ids_path to resolve variant_id.")
      variant_ids <- .load_variant_ids(variant_ids_path, ids_key)
    }
    idx <- match(variant_id, as.character(variant_ids))
    if (is.na(idx)) stop("variant_id '", variant_id, "' not found.")
  }
  idx <- as.integer(idx); if (idx < 1L) stop("idx must be >=1")
  
  prof1_key <- paste0("/", group, "/allele1_pred_profiles")
  prof2_key <- paste0("/", group, "/allele2_pred_profiles")
  cnt1_key  <- paste0("/", group, "/allele1_pred_counts")
  cnt2_key  <- paste0("/", group, "/allele2_pred_counts")
  
  prof1 <- .h5_read_profile_1variant(h5_file_path, prof1_key, idx)  # REF by convention
  prof2 <- .h5_read_profile_1variant(h5_file_path, prof2_key, idx)  # ALT by convention
  L <- length(prof1)
  if (L != length(prof2)) stop("Profile length mismatch: ", L, " vs ", length(prof2))
  
  cnt1 <- .h5_read_count_scalar(h5_file_path, cnt1_key, idx)       # REF total (or scalar)
  cnt2 <- .h5_read_count_scalar(h5_file_path, cnt2_key, idx)       # ALT total (or scalar)
  
  list(
    index = idx,
    L = L,
    ref_profile = prof1,
    alt_profile = prof2,
    ref_count = as.numeric(cnt1),
    alt_count = as.numeric(cnt2),
    delta_profile = prof2 - prof1,
    delta_count = as.numeric(cnt2) - as.numeric(cnt1),
    log2_ratio_alt_ref = if (cnt1 > 0) log2(as.numeric(cnt2) / as.numeric(cnt1)) else NA_real_
  )
}

# ---------- Seqinfo/BigWig export (robust, no seqlevels mismatch) ----------
# Build a Seqinfo from a chrom.sizes file (2 columns: name, length)
seqinfo_from_chrom_sizes <- function(chrom_sizes_path) {
  cs <- fread(chrom_sizes_path, col.names = c("seqnames","seqlengths"))
  Seqinfo(seqnames = cs$seqnames, seqlengths = setNames(as.numeric(cs$seqlengths), cs$seqnames))
}

# Export a numeric profile vector as 1-bp steps to BigWig
export_profile_bigwig <- function(scores, chr, start_pos, chrom_sizes_path, out_path) {
  scores <- as.numeric(scores)
  
  # seqinfo
  cs <- fread(chrom_sizes_path, col.names = c("seqnames","seqlengths"))
  si <- Seqinfo(seqnames = cs$seqnames, seqlengths = setNames(as.numeric(cs$seqlengths), cs$seqnames))
  
  if (!(chr %in% seqlevels(si))) {
    if (paste0("chr", chr) %in% seqlevels(si)) chr <- paste0("chr", chr)
    else if (sub("^chr", "", chr) %in% seqlevels(si)) chr <- sub("^chr", "", chr)
    else if (chr == "MT" && "chrM" %in% seqlevels(si)) chr <- "chrM"
    else if (chr == "chrM" && "MT" %in% seqlevels(si)) chr <- "MT"
  }
  
  # 🔴 THIS is where errors usually happen:
  end_pos <- start_pos + length(scores) - 1L
  
  gr <- GRanges(
    seqnames = chr,
    ranges   = IRanges(start = start_pos, end = end_pos),
    score    = scores    # score must be same length as ranges!
  )
  
  seqinfo(gr) <- si[chr]
  export.bw(gr, out_path)
}




# ---------- High-level: read from H5 and write BigWigs ----------
# If you already know chr/pos, you can align the profile around the variant position.
# Otherwise, you can pass a window start explicitly.
make_pred_profile_bigwigs_from_h5 <- function(
    h5_file_path,
    idx = NULL,
    variant_id = NULL,
    variant_ids = NULL,
    variant_ids_path = NULL,
    ids_key = "variant_ids",
    group = "observed",
    chr,
    pos,                        # 1-based variant position in the genome
    chrom_sizes_path,
    out_dir = ".",
    prefix = "pred",
    export_delta = TRUE,
    center_idx = NULL          # if NULL, center at L/2; window_start = pos - (center_idx-1)
) {
  pars <- parse_variant_from_h5(h5_file_path, idx, variant_id, variant_ids, variant_ids_path, ids_key, group)
  L <- pars$L
  if (is.null(center_idx)) center_idx <- ceiling(L/2)
  window_start <- as.integer(pos) - (center_idx - 1L)
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  base <- sprintf("%s_%s_%d", prefix, sub("^chr","", chr), as.integer(pos))
  
  ref_path   <- file.path(out_dir, sprintf("%s.ref.bw",   base))
  alt_path   <- file.path(out_dir, sprintf("%s.alt.bw",   base))
  delta_path <- if (export_delta) file.path(out_dir, sprintf("%s.delta.bw", base)) else NULL
  
  pars$ref_profile   <- as.numeric(pars$ref_profile)
  pars$alt_profile   <- as.numeric(pars$alt_profile)
  pars$delta_profile <- if (export_delta) as.numeric(pars$delta_profile) else NULL
  
  
  export_profile_bigwig(pars$ref_profile, chr, window_start, chrom_sizes_path, ref_path)
  export_profile_bigwig(pars$alt_profile, chr, window_start, chrom_sizes_path, alt_path)
  if (export_delta) export_profile_bigwig(pars$delta_profile, chr, window_start, chrom_sizes_path, delta_path)
  
  invisible(list(
    index         = pars$index,
    L             = L,
    center_idx    = center_idx,
    window_start  = window_start,
    ref_path      = ref_path,
    alt_path      = alt_path,
    delta_path    = delta_path,
    ref_count     = pars$ref_count,
    alt_count     = pars$alt_count,
    log2_ratio    = pars$log2_ratio_alt_ref
  ))
}

