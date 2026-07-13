#function to extract the allele prediction and shap contribution scores for a specific variant
#this doest work very vell
# ---- 1) load reticulate and initialize the python helper once ----
Sys.setenv(RETICULATE_PYTHON = "~/.conda/envs/renviron_ne1/bin/python")
library(reticulate)
py_config()

.init_variant_bw_exporter <- function(force = TRUE) {
  try(reticulate::py_run_string("try: del extract_allele_shap\nexcept NameError: pass"), silent = TRUE)
  reticulate::py_run_string("import os
import numpy as np
import pandas as pd

# deps that must be installed in your Python env:
#   pip install pyBigWig deepdish h5py
import pyBigWig
import deepdish as dd
import h5py

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

    # trim left
    if start0 < 0:
        values = values[-start0:]
        start0 = 0

    # trim right
    end0 = start0 + values.shape[0]
    if end0 > chrom_len:
        values = values[:chrom_len - start0]
        end0 = start0 + values.shape[0]

    if values.size == 0:
        raise ValueError('Nothing to write: empty vector after trimming to chrom bounds')

    header = [(c, int(chrom_sizes[c])) for c in chrom_sizes]
    bw = pyBigWig.open(bw_path, 'w')
    bw.addHeader(header)
    starts = list(range(start0, end0))
    ends   = [s+1 for s in starts]
    bw.addEntries([chrom]*len(starts), starts, ends=ends, values=values.tolist())
    bw.close()
    return bw_path

def _detect_snp_index_from_raw(raw_ref, raw_alt):
    if raw_ref is None or raw_alt is None:
        return None
    raw_ref = np.asarray(raw_ref)
    raw_alt = np.asarray(raw_alt)
    assert raw_ref.ndim == 2 and raw_alt.ndim == 2 and raw_ref.shape == raw_alt.shape and raw_ref.shape[0] == 4
    L = raw_ref.shape[1]
    diff = np.abs(raw_ref - raw_alt).sum(axis=0)  # (L,)
    cand = np.where(diff > 0)[0]
    if cand.size == 0:
        return L // 2
    return int(cand[np.argmax(diff[cand])])

def export_variant_bigwigs(data, atrr, tsv_df, variant_id, chrom_sizes_path,
                           outdir='.', prefix=None, write_channels=False):
    os.makedirs(outdir, exist_ok=True)
    chrom_sizes = _load_chrom_sizes(chrom_sizes_path)

    # TSV lookup; pandas keeps original row index as .index
    row = tsv_df.loc[tsv_df['variant_id'] == variant_id]
    if row.empty:
        raise ValueError(f'variant_id {variant_id} not found in TSV')
    row = row.iloc[0]
    i = int(row.name)  # aligns to data['observed'] order

    chrom = str(row['chr'])
    pos_1b = int(row['pos'])

    # predictions
    obs = data['observed']
    A1 = np.asarray(obs['allele1_pred_profiles'])  # (N, L_pred)
    A2 = np.asarray(obs['allele2_pred_profiles'])  # (N, L_pred)
    pred_ref = A1[i].squeeze()
    pred_alt = A2[i].squeeze()
    L_pred   = int(pred_ref.shape[0])

    center_idx_pred = L_pred // 2
    start0_pred = (pos_1b - 1) - center_idx_pred

    # attributions
    v_ids   = np.asarray(atrr['variant_ids'])
    alleles = np.asarray(atrr['alleles'])             # 0=allele1(REF), 1=allele2(ALT)
    SHAP    = np.asarray(atrr['projected_shap']['seq'])  # (N_attr, 4, L_attr)

    mask = (v_ids == variant_id)
    idx_ref = np.where(mask & (alleles == 0))[0]
    idx_alt = np.where(mask & (alleles == 1))[0]
    if idx_ref.size == 0 or idx_alt.size == 0:
        raise ValueError(f'Missing REF/ALT attributions for {variant_id}')

    ref_attr_full = np.asarray(SHAP[idx_ref[0]], dtype=np.float64)  # (4, L_attr)
    alt_attr_full = np.asarray(SHAP[idx_alt[0]], dtype=np.float64)  # (4, L_attr)
    L_attr = int(min(ref_attr_full.shape[-1], alt_attr_full.shape[-1]))
    ref_attr_full = ref_attr_full[:, :L_attr]
    alt_attr_full = alt_attr_full[:, :L_attr]

    # optional raw to nail exact SNP column
    snp_idx_attr = None
    if 'raw' in atrr and isinstance(atrr['raw'], dict) and 'seq' in atrr['raw']:
        raw_seq = np.asarray(atrr['raw']['seq'])
        raw_ref = raw_seq[idx_ref[0]]
        raw_alt = raw_seq[idx_alt[0]]
        snp_idx_attr = _detect_snp_index_from_raw(raw_ref, raw_alt)
    center_idx_attr = int(snp_idx_attr) if snp_idx_attr is not None else (L_attr // 2)
    center_idx_attr = max(0, min(center_idx_attr, L_attr-1))
    start0_attr = (pos_1b - 1) - center_idx_attr

    # match prediction length to attribution (for consistent track lengths)
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

    contrib_ref_sum = ref_attr_full.sum(axis=0)  # (L_attr,)
    contrib_alt_sum = alt_attr_full.sum(axis=0)  # (L_attr,)

    base = prefix if prefix else variant_id
    bw_pred_ref = os.path.join(outdir, f'{base}.pred_REF.bw')
    bw_pred_alt = os.path.join(outdir, f'{base}.pred_ALT.bw')
    bw_con_ref  = os.path.join(outdir, f'{base}.contribSum_REF.bw')
    bw_con_alt  = os.path.join(outdir, f'{base}.contribSum_ALT.bw')

    _write_bigwig_from_vector(chrom, start0_pred, pred_ref, bw_pred_ref, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_pred, pred_alt, bw_pred_alt, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_attr, contrib_ref_sum, bw_con_ref, chrom_sizes)
    _write_bigwig_from_vector(chrom, start0_attr, contrib_alt_sum, bw_con_alt, chrom_sizes)

    bw_channels = {}
    if write_channels:
        bases = ['A','C','G','T']
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
    # load deepdish dicts and TSV
    data = dd.io.load(pred_h5)
    atrr = dd.io.load(attr_h5)
    tsv_df = pd.read_csv(tsv_path, sep='\\t', header=0)
    # keep original 0..N-1 row indices (no reset) so .name is the row number
    return export_variant_bigwigs(data, atrr, tsv_df, variant_id, chrom_sizes_path,
                                  outdir=outdir, prefix=prefix, write_channels=write_channels)
  ")
  invisible(TRUE)
}

# # ---- 2) R wrapper you can call directly ----
# export_variant_bigwigs_py <- function(variant_id,
#                                       tsv_path,
#                                       pred_h5,
#                                       attr_h5,
#                                       chrom_sizes,   # path to UCSC-style chrom.sizes
#                                       outdir = ".",
#                                       prefix = NULL,
#                                       write_channels = FALSE) {
#   # ensure python code is loaded
# if (!py_exists("export_variant_bigwigs_from_paths")) {
#   init_variant_bw_exporter()
# }
# py$export_variant_bigwigs_from_paths(
#   pred_h5 = pred_h5,
#   attr_h5 = attr_h5,
#   tsv_path = tsv_path,
#   variant_id = variant_id,
#   chrom_sizes_path = chrom_sizes,
#   outdir = outdir,
#   prefix = if (is.null(prefix)) NULL else prefix,
#   write_channels = write_channels
# )
# }
