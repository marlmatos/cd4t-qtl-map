# --- deps ---
library(rhdf5)

# Flexible scalar read for counts datasets with shapes (1,N), (N,1), or (N)
.h5_read_count_scalar <- function(h5_path, dset, idx) {
  # Try (1, idx)
  out <- try(h5read(h5_path, dset, index = list(1L, as.integer(idx))), silent = TRUE)
  if (!inherits(out, "try-error")) return(as.numeric(out))
  # Try (idx, 1)
  out <- try(h5read(h5_path, dset, index = list(as.integer(idx), 1L)), silent = TRUE)
  if (!inherits(out, "try-error")) return(as.numeric(out))
  # Try 1D
  out <- try(h5read(h5_path, dset, index = list(as.integer(idx))), silent = TRUE)
  if (!inherits(out, "try-error")) return(as.numeric(out))
  stop(sprintf("Unable to read scalar from '%s' at index %s; dataset shape not supported.", dset, idx))
}

# Load variant_ids from a helper source (H5 or text)
.load_variant_ids <- function(src_path, ids_key = "variant_ids") {
  if (grepl("\\.h5$|\\.hdf5$", src_path, ignore.case = TRUE)) {
    if (!H5Lexists(src_path, ids_key)) {
      # Try with leading slash
      alt <- if (startsWith(ids_key, "/")) substring(ids_key, 2L) else paste0("/", ids_key)
      if (!H5Lexists(src_path, alt))
        stop(sprintf("Dataset '%s' not found in %s", ids_key, src_path))
      ids_key <- alt
    }
    ids <- h5read(src_path, ids_key)
    return(as.character(ids))
  } else {
    # Assume plain text, one ID per line (or TSV first column)
    con <- file(src_path, "r"); on.exit(close(con), add = TRUE)
    lines <- readLines(con)
    ids <- sub("\t.*$", "", lines)  # keep first column if TSV
    ids <- ids[nzchar(ids)]
    return(ids)
  }
}

# ---- Public helper ----
# Extract allelic prediction *counts* for a single variant from an averaged predictions H5.
extract_allelic_pred_counts_R <- function(
  h5_file_path,
  target_variant_id,
  variant_ids = NULL,        # character vector of IDs in the averaged file's order
  variant_ids_path = NULL,   # path to H5 (with dataset) or text file listing IDs
  ids_key = "variant_ids",   # dataset name in the IDs H5, if variant_ids_path is H5
  group = "observed"         # group name in averaged file
) {
  stopifnot(file.exists(h5_file_path))

  # Resolve IDs
  if (is.null(variant_ids)) {
    if (is.null(variant_ids_path)) {
      stop("Averaged file typically lacks variant_ids. Provide either `variant_ids` vector ",
           "or `variant_ids_path` (H5 dataset or text file).")
    }
    variant_ids <- .load_variant_ids(variant_ids_path, ids_key)
  }
  variant_ids <- as.character(variant_ids)

  # Find index for the target
  idx <- match(target_variant_id, variant_ids)
  if (is.na(idx)) {
    stop(sprintf("target_variant_id '%s' not found in provided variant_ids.", target_variant_id))
  }

  # Dataset paths
  d1 <- paste0("/", group, "/allele1_pred_counts")
  d2 <- paste0("/", group, "/allele2_pred_counts")

  # Read scalars robustly
  a1 <- .h5_read_count_scalar(h5_file_path, d1, idx)
  a2 <- .h5_read_count_scalar(h5_file_path, d2, idx)

  # Return a clean list (include convenience deltas)
  list(
    variant_id = target_variant_id,
    index      = as.integer(idx),
    allele1_pred_count = a1,
    allele2_pred_count = a2,
    delta_alt_minus_ref = a2 - a1,
    log2_ratio_alt_ref  = if (a1 > 0) log2(a2 / a1) else NA_real_
  )
}

# ------------- Example usage -------------
# averaged_h5 <- "/gpfs/.../averaged_variant_prediction_scores.h5"
# ids_txt_or_h5 <- "/gpfs/.../variant_ids.txt"  # or an H5 with dataset 'variant_ids'
# res <- extract_allelic_pred_counts_R(
#   h5_file_path = averaged_h5,
#   target_variant_id = "chr1_153617778_C_A",
#   variant_ids_path = ids_txt_or_h5,   # or pass `variant_ids = your_vector`
#   ids_key = "variant_ids",
#   group = "observed"
# )
# str(res)
