# ============================================================================
# ChromBPNet Variant BigWig Export Helper
# Complete version with initialization and safe wrappers
# ============================================================================
# ---- 1) load reticulate and initialize the python helper once ----
Sys.setenv(RETICULATE_PYTHON = "~/.conda/envs/renviron_ne1/bin/python")
library(reticulate)
py_config()

# Load reticulate if not already loaded
if (!require("reticulate", quietly = TRUE)) {
  stop("Package 'reticulate' is required. Install with: install.packages('reticulate')")
}

# ============================================================================
# PART 1: Python Function Initialization
# ============================================================================

.init_variant_bw_exporter <- function(force = TRUE) {
  # Force delete ALL related Python functions to avoid caching issues
  try(py_run_string("
try: 
    del export_variant_bigwigs_from_paths
    del export_variant_bigwigs
    del _load_attr_with_variant_ids
    del _load_chrom_sizes
    del _write_bigwig_from_vector
    del _detect_snp_index_from_raw
except NameError: 
    pass
"), silent = TRUE)
  
  # Now define everything fresh
  py_run_string("
import os
import numpy as np
import pandas as pd
import pyBigWig
import deepdish as dd
import h5py
import pickle

def _load_chrom_sizes(chrom_sizes_path):
    sizes = {}
    with open(chrom_sizes_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, sz = line.rstrip().split('\\t')[:2]
            sizes[chrom] = int(sz)
    return sizes

def _write_bigwig_from_vector(chrom, start0, values, bw_path, chrom_sizes):
    if chrom not in chrom_sizes:
        raise ValueError(f'Chrom {chrom} not in chrom sizes')
    chrom_len = chrom_sizes[chrom]
    start0 = int(start0)
    values = np.asarray(values, dtype=float)
    if start0 < 0:
        values = values[-start0:]
        start0 = 0
    end0 = start0 + values.shape[0]
    if end0 > chrom_len:
        values = values[:chrom_len - start0]
        end0 = start0 + values.shape[0]
    if values.size == 0:
        raise ValueError('Nothing to write: empty vector after trimming')
    header = [(c, int(chrom_sizes[c])) for c in chrom_sizes]
    bw = pyBigWig.open(bw_path, 'w')
    bw.addHeader(header)
    starts = list(range(start0, end0))
    ends = [s + 1 for s in starts]
    bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values.tolist())
    bw.close()
    return bw_path

def _detect_snp_index_from_raw(raw_ref, raw_alt):
    if raw_ref is None or raw_alt is None:
        return None
    raw_ref = np.asarray(raw_ref)
    raw_alt = np.asarray(raw_alt)
    assert raw_ref.ndim == 2 and raw_alt.ndim == 2
    L = raw_ref.shape[1]
    diff = np.abs(raw_ref - raw_alt).sum(axis=0)
    cand = np.where(diff > 0)[0]
    if cand.size == 0:
        return L // 2
    return int(cand[np.argmax(diff[cand])]) 

def _load_attr_with_variant_ids(attr_h5):
    '''Load attribution H5 and unpickle variant_ids.'''
    atrr = dd.io.load(attr_h5)
    
    # Always reload variant_ids via pickle
    with h5py.File(attr_h5, 'r') as f:
        if 'variant_ids' not in f:
            raise KeyError('variant_ids dataset not found in H5')
        node = f['variant_ids']
        raw = node[0]
        blob = bytes(raw)
        variant_ids = pickle.loads(blob)
    
    atrr['variant_ids'] = variant_ids
    return atrr

def export_variant_bigwigs(data, atrr, tsv_df, variant_id, chrom_sizes_path,
                           outdir='.', prefix=None, write_channels=False):
    os.makedirs(outdir, exist_ok=True)
    chrom_sizes = _load_chrom_sizes(chrom_sizes_path)
    
    # TSV lookup
    row = tsv_df.loc[tsv_df['variant_id'] == variant_id]
    if row.empty:
        raise ValueError(f'variant_id {variant_id} not found in TSV')
    row = row.iloc[0]
    i = int(row.name)
    chrom = str(row['chr'])
    pos_1b = int(row['pos'])
    
    # Predictions
    obs = data['observed']
    A1 = np.asarray(obs['allele1_pred_profiles'])
    A2 = np.asarray(obs['allele2_pred_profiles'])
    pred_ref = A1[i].squeeze()
    pred_alt = A2[i].squeeze()
    L_pred = int(pred_ref.shape[0])
    center_idx_pred = L_pred // 2
    start0_pred = (pos_1b - 1) - center_idx_pred
    
    # Attributions
    variant_id = str(variant_id)
    v_ids = np.asarray(atrr['variant_ids'])
    
    # Normalize to strings
    if v_ids.dtype.kind in ('S', 'U', 'O'):
        v_ids = v_ids.astype(str)
    else:
        v_ids = np.vectorize(str)(v_ids)
    
    alleles = np.asarray(atrr['alleles'])
    SHAP = np.asarray(atrr['projected_shap']['seq'])
    
    mask = (v_ids == variant_id)
    
    if not np.any(mask):
        example = list(v_ids[:5])
        raise ValueError(
            f'variant_id {variant_id} not found in attribution variant_ids. '
            f'Example IDs (first 5): {example}'
        )
    
    idx_ref = np.where(mask & (alleles == 0))[0]
    idx_alt = np.where(mask & (alleles == 1))[0]
    
    if idx_ref.size == 0 or idx_alt.size == 0:
        raise ValueError(
            f'Found variant_id {variant_id} but missing REF/ALT alleles. '
            f'alleles present: {alleles[mask]}'
        )
    
    ref_attr_full = np.asarray(SHAP[idx_ref[0]], dtype=np.float64)
    alt_attr_full = np.asarray(SHAP[idx_alt[0]], dtype=np.float64)
    L_attr = int(min(ref_attr_full.shape[-1], alt_attr_full.shape[-1]))
    ref_attr_full = ref_attr_full[:, :L_attr]
    alt_attr_full = alt_attr_full[:, :L_attr]
    
    # SNP detection
    snp_idx_attr = None
    if 'raw' in atrr and isinstance(atrr['raw'], dict) and 'seq' in atrr['raw']:
        raw_seq = np.asarray(atrr['raw']['seq'])
        raw_ref = raw_seq[idx_ref[0]]
        raw_alt = raw_seq[idx_alt[0]]
        snp_idx_attr = _detect_snp_index_from_raw(raw_ref, raw_alt)
    center_idx_attr = int(snp_idx_attr) if snp_idx_attr is not None else (L_attr // 2)
    center_idx_attr = max(0, min(center_idx_attr, L_attr - 1))
    start0_attr = (pos_1b - 1) - center_idx_attr
    
    # Length matching
    if L_pred != L_attr:
        if L_attr <= L_pred:
            left = max(0, center_idx_pred - (L_attr // 2))
            right = left + L_attr
            pred_ref = pred_ref[left:right]
            pred_alt = pred_alt[left:right]
            center_idx_pred = L_attr // 2
            start0_pred = (pos_1b - 1) - center_idx_pred
        else:
            pad_left = (L_attr - L_pred) // 2
            pad_right = L_attr - L_pred - pad_left
            pred_ref = np.pad(pred_ref, (pad_left, pad_right), constant_values=0.0)
            pred_alt = np.pad(pred_alt, (pad_left, pad_right), constant_values=0.0)
            center_idx_pred = L_attr // 2
            start0_pred = (pos_1b - 1) - center_idx_pred
    
    contrib_ref_sum = ref_attr_full.sum(axis=0)
    contrib_alt_sum = alt_attr_full.sum(axis=0)
    
    base = prefix if prefix else str(variant_id)
    bw_pred_ref = os.path.join(outdir, f'{base}.pred_REF.bw')
    bw_pred_alt = os.path.join(outdir, f'{base}.pred_ALT.bw')
    bw_con_ref = os.path.join(outdir, f'{base}.contribSum_REF.bw')
    bw_con_alt = os.path.join(outdir, f'{base}.contribSum_ALT.bw')
    
    _write_bigwig_from_vector(chrom, start0_pred, pred_ref, bw_pred_ref, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_pred, pred_alt, bw_pred_alt, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_attr, contrib_ref_sum, bw_con_ref, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_attr, contrib_alt_sum, bw_con_alt, chrom_sizes)
    
    bw_channels = {}
    if write_channels:
        bases = ['A', 'C', 'G', 'T']
        for bi, b in enumerate(bases):
            p_ref = os.path.join(outdir, f'{base}.contrib_{b}_REF.bw')
            p_alt = os.path.join(outdir, f'{base}.contrib_{b}_ALT.bw')
            _write_bigwig_from_vector(chrom, start0_attr, ref_attr_full[bi, :], p_ref, chrom_sizes)
            _write_bigwig_from_vector(chrom, start0_attr, alt_attr_full[bi, :], p_alt, chrom_sizes)
            bw_channels[b + '_REF'] = p_ref
            bw_channels[b + '_ALT'] = p_alt
    
    return {
        'pred_REF': bw_pred_ref,
        'pred_ALT': bw_pred_alt,
        'contribSum_REF': bw_con_ref,
        'contribSum_ALT': bw_con_alt,
        'channels': bw_channels,
        'centers': {'pred_center_idx': int(center_idx_pred), 'attr_center_idx': int(center_idx_attr)}
    }

def export_variant_bigwigs_from_paths(pred_h5, attr_h5, tsv_path, variant_id,
                                      chrom_sizes_path, outdir='.', prefix=None, write_channels=False):
    data = dd.io.load(pred_h5)
    atrr = _load_attr_with_variant_ids(attr_h5)
    tsv_df = pd.read_csv(tsv_path, sep='\\t', header=0)
    return export_variant_bigwigs(data, atrr, tsv_df, variant_id, chrom_sizes_path,
                                  outdir=outdir, prefix=prefix, write_channels=write_channels)
")
  
  message("✓ Python functions initialized successfully")
  invisible(TRUE)
}

# ============================================================================
# PART 2: Safe R Wrapper Functions
# ============================================================================

export_variant_bigwigs_safe <- function(
    pred_h5,
    attr_h5, 
    tsv_path,
    variant_id,
    chrom_sizes_path,
    outdir = NULL,
    prefix = NULL,
    write_channels = FALSE
) {
  
  # Check if Python functions are loaded
  if (!py_has_attr(py, "export_variant_bigwigs_from_paths")) {
    stop(
      "Python function 'export_variant_bigwigs_from_paths' not found.\n",
      "You must run .init_variant_bw_exporter() first!"
    )
  }
  
  # Input validation
  if (!file.exists(pred_h5)) {
    stop("pred_h5 file not found: ", pred_h5)
  }
  if (!file.exists(attr_h5)) {
    stop("attr_h5 file not found: ", attr_h5)
  }
  if (!file.exists(tsv_path)) {
    stop("tsv_path file not found: ", tsv_path)
  }
  if (!file.exists(chrom_sizes_path)) {
    stop("chrom_sizes_path file not found: ", chrom_sizes_path)
  }
  
  # Handle output directory
  if (is.null(outdir) || outdir == "") {
    outdir <- getwd()
    message("No outdir specified, using current directory: ", outdir)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    message("Creating output directory: ", outdir)
    dir.create(outdir, recursive = TRUE)
  }
  
  # Try the export with error handling
  result <- tryCatch(
    {
      py$export_variant_bigwigs_from_paths(
        pred_h5 = pred_h5,
        attr_h5 = attr_h5,
        tsv_path = tsv_path,
        variant_id = variant_id,
        chrom_sizes_path = chrom_sizes_path,
        outdir = outdir,
        prefix = prefix,
        write_channels = write_channels
      )
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      
      # Check for common errors
      if (grepl("not found in TSV", err_msg)) {
        stop(sprintf(
          "Variant '%s' not found in TSV file.\nCheck variant ID format in: %s",
          variant_id, tsv_path
        ))
      } else if (grepl("not found in attribution", err_msg)) {
        stop(sprintf(
          "Variant '%s' found in TSV but not in attribution file.\nThis may indicate mismatched files: %s",
          variant_id, attr_h5
        ))
      } else if (grepl("missing REF/ALT alleles", err_msg)) {
        stop(sprintf(
          "Variant '%s' found but missing REF or ALT allele in attribution file.\nThis indicates incomplete attribution data for this variant.",
          variant_id
        ))
      } else {
        # Re-throw other errors
        stop(sprintf("Error processing variant '%s':\n%s", variant_id, err_msg))
      }
    }
  )
  
  # Success message
  if (!is.null(result)) {
    message(sprintf("\n✓ Successfully exported bigWigs for variant: %s", variant_id))
    message(sprintf("  Output directory: %s", outdir))
    if (!is.null(result$pred_REF)) {
      message(sprintf("  Files created:"))
      message(sprintf("    - %s", basename(result$pred_REF)))
      message(sprintf("    - %s", basename(result$pred_ALT)))
      message(sprintf("    - %s", basename(result$contribSum_REF)))
      message(sprintf("    - %s", basename(result$contribSum_ALT)))
    }
  }
  
  invisible(result)
}

# Batch export function
export_variants_batch <- function(
    pred_h5,
    attr_h5,
    tsv_path,
    variant_ids,
    chrom_sizes_path,
    outdir = NULL,
    write_channels = FALSE,
    skip_on_error = TRUE
) {
  
  if (is.null(outdir) || outdir == "") {
    outdir <- getwd()
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  n_variants <- length(variant_ids)
  message(sprintf("Processing %d variants...\n", n_variants))
  
  results <- list()
  failed <- character()
  
  for (i in seq_along(variant_ids)) {
    vid <- variant_ids[i]
    message(sprintf("[%d/%d] Processing: %s", i, n_variants, vid))
    
    result <- tryCatch(
      {
        export_variant_bigwigs_safe(
          pred_h5 = pred_h5,
          attr_h5 = attr_h5,
          tsv_path = tsv_path,
          variant_id = vid,
          chrom_sizes_path = chrom_sizes_path,
          outdir = outdir,
          prefix = NULL,
          write_channels = write_channels
        )
      },
      error = function(e) {
        if (skip_on_error) {
          warning(sprintf("Failed to process %s: %s", vid, conditionMessage(e)))
          failed <<- c(failed, vid)
          NULL
        } else {
          stop(e)
        }
      }
    )
    
    if (!is.null(result)) {
      results[[vid]] <- result
    }
  }
  
  # Summary
  message(sprintf("\n=== Batch Export Summary ==="))
  message(sprintf("Total variants: %d", n_variants))
  message(sprintf("Successful: %d", length(results)))
  message(sprintf("Failed: %d", length(failed)))
  
  if (length(failed) > 0) {
    message("\nFailed variants:")
    for (f in failed) {
      message(sprintf("  - %s", f))
    }
  }
  
  invisible(list(
    results = results,
    failed = failed,
    success_rate = length(results) / n_variants
  ))
}

# ============================================================================
# Usage message
# ============================================================================

cat("
================================================================================
ChromBPNet Variant BigWig Exporter - COMPLETE VERSION
================================================================================

STEP 1: Initialize Python functions (do this once per session)
---------------------------------------------------------------
.init_variant_bw_exporter()


STEP 2: Export variants
------------------------

Single variant:
~~~~~~~~~~~~~~~
result <- export_variant_bigwigs_safe(
  pred_h5 = '/path/to/predictions.h5',
  attr_h5 = '/path/to/attributions.h5',
  tsv_path = '/path/to/variants.tsv',
  variant_id = 'chr6_29708854_G_A',
  chrom_sizes_path = '/path/to/chrom.sizes',
  outdir = '/path/to/output',
  write_channels = FALSE
)

Batch export:
~~~~~~~~~~~~~
results <- export_variants_batch(
  pred_h5 = '/path/to/predictions.h5',
  attr_h5 = '/path/to/attributions.h5',
  tsv_path = '/path/to/variants.tsv',
  variant_ids = c('chr6_29708854_G_A', 'chr1_1000112_G_T'),
  chrom_sizes_path = '/path/to/chrom.sizes',
  outdir = '/path/to/output',
  skip_on_error = TRUE
)

================================================================================
")