# --- deps ---
library(reticulate)

.init_py_extract_allele_shap <- function(force = TRUE) {
  # nuke any old def so we don’t keep stale code around
  try(reticulate::py_run_string("try: del extract_allele_shap\nexcept NameError: pass"), silent = TRUE)
  
  reticulate::py_run_string("
import numpy as np
try:
    import deepdish as dd
except Exception as e:
    raise ImportError('Need deepdish: pip install deepdish') from e

def _find_axis_of_size(shape, size):
    return [i for i, s in enumerate(shape) if s == size]

def _pick_axis(shape, candidates, name):
    if not candidates:
        raise ValueError(f'Could not find {name} axis in tensor of shape {shape}')
    if len(candidates) > 1:
        return candidates[-1]
    return candidates[0]

def extract_allele_shap(h5_file_path, target_variant_id,
                        variant_ids_key='variant_ids', shap_key='shap', subkey='seq'):
    data = dd.io.load(h5_file_path)

    if variant_ids_key not in data:
        raise KeyError(f'Missing key `{variant_ids_key}` in H5')
    variant_ids = data[variant_ids_key]
    try:
        variant_ids = variant_ids.astype(str)
    except Exception:
        variant_ids = np.array([str(v) for v in variant_ids])
    n_vars = int(len(variant_ids))

    shap = data[shap_key][subkey]
    shape = tuple(shap.shape)

    allele_axes = _find_axis_of_size(shape, 4)
    var_axes    = _find_axis_of_size(shape, n_vars)
    allele_ax = _pick_axis(shape, allele_axes, 'allele(4)')
    var_ax    = _pick_axis(shape, var_axes, 'variant(n_vars)')

    remaining = [i for i in range(len(shape)) if i not in (allele_ax, var_ax)]
    if not remaining:
        raise ValueError(f'Cannot infer length axis from shape {shape}')
    cand_len = [i for i in remaining if shape[i] >= 1000] or remaining
    len_ax   = cand_len[0]
    seq_len  = int(shape[len_ax])

    matches = np.where(variant_ids == str(target_variant_id))[0]
    if matches.size == 0:
        raise ValueError(f'Variant ID not found: {target_variant_id}')
    vidx = int(matches[0])

    parts = str(target_variant_id).split('_')
    if len(parts) != 4:
        raise ValueError('Expected variant_id format chr_pos_ref_alt')
    chr_name, pos_s, ref_allele, alt_allele = parts
    try:
        pos = int(pos_s)
    except Exception:
        pos = pos_s
    base_to_idx = {'A':0, 'C':1, 'G':2, 'T':3}
    if ref_allele not in base_to_idx or alt_allele not in base_to_idx:
        raise ValueError('Alleles must be A/C/G/T')

    indexer = [slice(None)] * len(shape)
    indexer[var_ax] = vidx

    def take_allele(base):
        ax = base_to_idx[base]
        idxr = list(indexer)
        idxr[allele_ax] = ax
        arr = shap[tuple(idxr)]
        # move length axis to last, reduce any remaining dims
        if arr.ndim > 1:
            arr = np.moveaxis(arr, len_ax if len_ax < len(shape) else -1, -1)
            if arr.ndim > 1:
                arr = np.reshape(arr, (-1, arr.shape[-1])).mean(axis=0)
        return np.asarray(arr, dtype=float)

    ref_shap = take_allele(ref_allele)
    alt_shap = take_allele(alt_allele)

    if ref_shap.shape[0] != seq_len or alt_shap.shape[0] != seq_len:
        raise RuntimeError(
            f'Length mismatch: seq_len={seq_len}, ref={ref_shap.shape}, alt={alt_shap.shape}; '
            f'shape={shape}, var_ax={var_ax}, allele_ax={allele_ax}, len_ax={len_ax}'
        )

    center_idx = int(seq_len // 2)

    return {
        'variant_id': str(target_variant_id),
        'chr': str(chr_name),
        'pos': pos,
        'ref_allele': str(ref_allele),
        'alt_allele': str(alt_allele),
        'ref_shap': ref_shap,
        'alt_shap': alt_shap,
        'seq_len': seq_len,
        'axes': {'allele_ax': int(allele_ax), 'var_ax': int(var_ax), 'len_ax': int(len_ax), 'shape': tuple(shape)},
        'center_idx': center_idx
    }
")
  invisible(TRUE)
}

extract_allele_shap_R <- function(h5_file_path, target_variant_id,
                                  variant_ids_key = "variant_ids",
                                  shap_key = "shap", subkey = "seq") {
  .init_py_extract_allele_shap()
  out <- reticulate::py$extract_allele_shap(h5_file_path, target_variant_id,
                                            variant_ids_key = variant_ids_key,
                                            shap_key = shap_key, subkey = subkey)
  out <- reticulate::py_to_r(out)
  # safe scalars
  pick1_chr <- function(x) as.character(x)[1]
  pick1_int <- function(x) as.integer(x[1])
  pick_num  <- function(x) as.numeric(x)
  
  list(
    variant_id = pick1_chr(out$variant_id),
    chr        = pick1_chr(out$chr),
    pos        = pick1_int(out$pos),
    ref_allele = pick1_chr(out$ref_allele),
    alt_allele = pick1_chr(out$alt_allele),
    seq_len    = pick1_int(out$seq_len),
    center_idx = pick1_int(out$center_idx),
    ref_shap   = pick_num(out$ref_shap),
    alt_shap   = pick_num(out$alt_shap),
    axes       = out$axes
  )
}


# --- example usage ---
# h5_file_path <- "/gpfs/commons/home/mmatos/cd4_chrombpnet/chrombpnet_model_b7/variant_contribution_scores/averaged_cd4_tcells_AJ_common_variants.shap.counts.h5"
# res <- extract_allele_shap_R(h5_file_path, "chr1_153617778_C_A")
# str(res)
# plot(res$ref_shap, type='l'); lines(res$alt_shap)



########################################
#######################################
########################################
make_shap_bigwig <- function(
    shap_data,
    ref_shap,
    alt_shap,
    chrom_sizes_path,
    out_dir = ".",
    prefix  = "shap",
    export_delta = TRUE,
    chr_style = c("UCSC","asis"),   # if your chrom.sizes use 'chr1' etc.
    center_idx = NULL               # if NULL, uses ceiling(L/2)
) {
  stopifnot(is.list(shap_data))
  chr <- shap_data$chr
  pos <- as.integer(shap_data$pos)
  ref <- as.character(shap_data$ref_allele)
  alt <- as.character(shap_data$alt_allele)
  
  ref_scores <- as.numeric(ref_shap)
  alt_scores <- as.numeric(alt_shap)
  stopifnot(length(ref_scores) == length(alt_scores))
  L <- length(ref_scores)
  if (is.null(center_idx)) center_idx <- ceiling(L/2)
  
  # --- chrom sizes & seqinfo ---
  chrom_sizes <- data.table::fread(chrom_sizes_path, col.names = c("chrom","size"))
  seqlens <- setNames(as.numeric(chrom_sizes$size), chrom_sizes$chrom)
  if (match.arg(chr_style) == "UCSC") {
    # ensure chr has 'chr' prefix like your sizes file
    if (!startsWith(chr, "chr") && any(startsWith(names(seqlens), "chr"))) {
      chr <- paste0("chr", chr)
    }
  }
  
  chr_len <- seqlens[[chr]]
  if (is.na(chr_len)) stop(sprintf("Chromosome %s not found in chrom.sizes", chr))
  
  # --- align SHAP window so center_idx sits at genomic position 'pos' ---
  start0 <- pos - (center_idx - 1L)
  end0   <- start0 + L - 1L
  
  # --- clip to chromosome bounds ---
  clip_start <- max(1L, start0)
  clip_end   <- min(chr_len, end0)
  
  left_clip  <- clip_start - start0
  right_clip <- end0 - clip_end
  idx_start  <- 1L + max(0L, left_clip)
  idx_end    <- L - max(0L, right_clip)
  
  if (idx_start > idx_end) stop("After clipping, no bases remain to export.")
  
  ref_scores_clip <- ref_scores[idx_start:idx_end]
  alt_scores_clip <- alt_scores[idx_start:idx_end]
  pos_vec <- seq.int(clip_start, clip_end)
  
  # --- helper to build GRanges ---
  make_gr <- function(scores, track_name_suffix) {
    gr <- GenomicRanges::GRanges(
      seqnames = chr,
      ranges   = IRanges::IRanges(start = pos_vec, width = 1L),
      score    = scores
    )
    # Set the seqlength ONLY for the seqlevel in `gr`
    GenomeInfoDb::seqlengths(gr) <- setNames(seqlens[chr], chr)  # << change here
    S4Vectors::mcols(gr)$name <- sprintf("%s_%s_%s_%s_%s",
                                         chr, pos, ref, alt, track_name_suffix)
    gr
  }
  
  gr_ref <- make_gr(ref_scores_clip, paste0("REF_", ref))
  gr_alt <- make_gr(alt_scores_clip, paste0("ALT_", alt))
  if (export_delta) {
    gr_dlt <- make_gr(alt_scores_clip - ref_scores_clip, "DELTA_ALTminusREF")
  }
  
  # --- ensure output dir exists ---
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # --- export BigWigs ---
  base <- sprintf("%s_%s_%d_%s_%s", prefix, sub("^chr","", chr), pos, ref, alt)
  out_ref_bw <- file.path(out_dir, sprintf("%s.ref.bw",  base))
  out_alt_bw <- file.path(out_dir, sprintf("%s.alt.bw",  base))
  rtracklayer::export.bw(gr_ref, out_ref_bw)
  rtracklayer::export.bw(gr_alt, out_alt_bw)
  
  out_dlt_bw <- NULL
  if (export_delta) {
    out_dlt_bw <- file.path(out_dir, sprintf("%s.delta.bw", base))
    rtracklayer::export.bw(gr_dlt, out_dlt_bw)
  }
  
  message(sprintf("Exported:\n  %s\n  %s%s",
                  out_ref_bw, out_alt_bw,
                  if (export_delta) paste0("\n  ", out_dlt_bw) else ""))
  
  invisible(list(
    ref_bw = out_ref_bw,
    alt_bw = out_alt_bw,
    delta_bw = out_dlt_bw,
    gr_ref = gr_ref,
    gr_alt = gr_alt,
    gr_delta = if (export_delta) gr_dlt else NULL,
    window = list(
      shap_len = L, center_idx = center_idx,
      start = clip_start, end = clip_end, pos = pos
    )
  ))
}


###usage

# res_bw <- make_shap_bigwig(
#   shap_data = shap_data,
#   ref_shap  = shap_data$ref_shap,
#   alt_shap  = shap_data$alt_shap,
#   chrom_sizes_path = "/gpfs/commons/home/mmatos/resources/genome/hg38.chrom.sizes",
#   out_dir   = "bw_out",
#   prefix    = "shap",
#   export_delta = TRUE
# )