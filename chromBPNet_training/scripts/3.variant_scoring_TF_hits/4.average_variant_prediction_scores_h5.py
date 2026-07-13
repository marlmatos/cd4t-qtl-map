#!/usr/bin/env python3

import os, sys, glob, json, argparse
import h5py
import numpy as np
import pandas as pd
import hdf5plugin


def _load_chrom_sizes(chrom_sizes_path):
    sizes = {}
    with open(chrom_sizes_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, sz = line.rstrip().split('\t')[:2]
            sizes[chrom] = int(sz)
    return sizes


def _read_all_datasets(h5, keep_keys=None, dataset_prefix=None):
    out = {}
    def visitor(name, obj):
        if not isinstance(obj, h5py.Dataset):
            return
        if keep_keys is not None:
            if name in keep_keys:
                out[name] = obj[()]
            return
        if dataset_prefix:
            if name.startswith(dataset_prefix + "/"):
                out[name] = obj[()]
            return
        out[name] = obj[()]
    h5.visititems(visitor)
    return out


def _group_by_key(list_of_dicts):
    out = {}
    for d in list_of_dicts:
        for k, v in d.items():
            out.setdefault(k, []).append(v)
    return out


def concat_chroms_within_fold(per_chrom_h5_paths, tsv_path,
                               keep_keys=None, dataset_prefix=None):
    """
    Load per-chrom H5s for one fold, attach variant_ids from TSV
    by matching chromosome, then concatenate along axis 0 (variants).
    """
    tsv = pd.read_csv(tsv_path, sep='\t', header=None,
                  names=['chr', 'pos', 'allele1', 'allele2', 'variant_id'])

    per_file = []
    variant_ids_per_file = []

    for p in sorted(per_chrom_h5_paths):
        # extract chrom from filename
        basename = os.path.basename(p)
        # expects something like ...chr1.variant_predictions.h5
        chrom = None
        for part in basename.replace('.', '_').split('_'):
            if part.startswith('chr'):
                chrom = part
                break
        if chrom is None:
            raise ValueError(f"Could not extract chrom from filename: {basename}")

        chrom_tsv = tsv[tsv['chr'] == chrom].reset_index(drop=True)
        if chrom_tsv.empty:
            raise ValueError(f"No variants in TSV for chrom {chrom}")

        with h5py.File(p, 'r') as f:
            data = _read_all_datasets(f, keep_keys=keep_keys,
                                       dataset_prefix=dataset_prefix)

        # validate row count matches
        n_variants = len(chrom_tsv)
        for k, arr in data.items():
            if arr.shape[0] != n_variants:
                raise ValueError(
                    f"{p}: dataset '{k}' has {arr.shape[0]} rows "
                    f"but TSV has {n_variants} variants for {chrom}"
                )

        per_file.append(data)
        variant_ids_per_file.append(chrom_tsv['variant_id'].values)

    # concatenate across chroms
    by_key = _group_by_key(per_file)
    fold_cat = {}
    for key, arrs in by_key.items():
        fold_cat[key] = np.concatenate(arrs, axis=0)

    fold_cat['variant_ids'] = np.concatenate(variant_ids_per_file, axis=0)
    return fold_cat


def align_folds_by_ids(concatenated_per_fold):
    """
    Align all folds to the variant order of fold 0, on axis 0.
    """
    ref_ids = concatenated_per_fold[0]['variant_ids'].astype(str)
    ref_map = {v: i for i, v in enumerate(ref_ids)}

    aligned = []
    for fi, d in enumerate(concatenated_per_fold):
        ids = d['variant_ids'].astype(str)

        if set(ids) != set(ref_ids):
            raise ValueError(
                f"Fold {fi} variant IDs do not match fold 0. "
                f"ref={len(ref_ids)} fold={len(ids)}"
            )

        sort_order = np.array([ref_map[v] for v in ids], dtype=int)

        d_aligned = {}
        for k, arr in d.items():
            if k == 'variant_ids':
                d_aligned[k] = arr
            elif arr.ndim >= 1 and arr.shape[0] == len(ids):
                d_aligned[k] = arr[sort_order]
            else:
                d_aligned[k] = arr
        aligned.append(d_aligned)

    return aligned


def average_folds(concatenated_per_fold, cast_back_float16=True):
    # exclude variant_ids from averaging
    keys = set(concatenated_per_fold[0].keys()) - {'variant_ids'}

    for d in concatenated_per_fold[1:]:
        missing = keys.symmetric_difference(set(d.keys()) - {'variant_ids'})
        if missing:
            raise ValueError(f"Key mismatch across folds: {missing}")

    averages = {}
    for key in sorted(keys):
        stacks = [d[key] for d in concatenated_per_fold]
        shape0 = stacks[0].shape
        if any(s.shape != shape0 for s in stacks[1:]):
            raise ValueError(
                f"Shape mismatch for key {key}: {[s.shape for s in stacks]}"
            )
        dtype0 = stacks[0].dtype
        acc = np.zeros(shape0,
                       dtype=(np.float32 if dtype0 == np.float16 else dtype0))
        for s in stacks:
            acc += s.astype(acc.dtype, copy=False)
        avg = acc / float(len(stacks))
        if cast_back_float16 and dtype0 == np.float16:
            avg = avg.astype(np.float16)
        averages[key] = avg

    # carry variant_ids from fold 0
    averages['variant_ids'] = concatenated_per_fold[0]['variant_ids']
    return averages


def save_h5_dict(h5_out_path, data_dict):
    os.makedirs(os.path.dirname(h5_out_path), exist_ok=True)
    with h5py.File(h5_out_path, 'w') as f:
        for key, arr in data_dict.items():
            if '/' in key:
                group, dset = key.split('/', 1)
            else:
                group, dset = '', key
            grp = f if group == '' else f.require_group(group)
            grp.create_dataset(dset, data=arr)
    return h5_out_path


def concat_then_average_folds(
    fold_glob,
    tsv_path,
    per_fold_file_glob="**/*.h5",
    out_path="averaged_folds.h5",
    keep_keys=None,
    dataset_prefix=None,
    dry_run=False,
):
    fold_dirs = sorted(glob.glob(fold_glob))
    if not fold_dirs:
        raise FileNotFoundError(f"No fold dirs matched: {fold_glob}")

    per_fold_concat = []
    for fd in fold_dirs:
        h5s = sorted(glob.glob(
            os.path.join(fd, per_fold_file_glob), recursive=True))
        if not h5s:
            print(f"[WARN] No H5 files in {fd}")
            continue
        print(f"[INFO] Fold {fd}: {len(h5s)} per-chrom files")
        fold_cat = concat_chroms_within_fold(
            h5s,
            tsv_path=tsv_path,
            keep_keys=keep_keys,
            dataset_prefix=dataset_prefix,
        )
        print(f"       total variants: {len(fold_cat['variant_ids'])}")
        per_fold_concat.append(fold_cat)

    if not per_fold_concat:
        raise RuntimeError("No per-fold data assembled.")

    print("[INFO] Aligning folds by variant_ids")
    per_fold_concat = align_folds_by_ids(per_fold_concat)

    print("[INFO] Averaging across folds")
    averaged = average_folds(per_fold_concat)

    if dry_run:
        print("[DRY-RUN] Would write:", out_path)
        for k, v in averaged.items():
            print(f"  {k}: shape={v.shape} dtype={v.dtype}")
        return out_path

    save_h5_dict(out_path, averaged)
    print(f"[OK] Saved averaged H5 to {out_path}")
    return out_path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Concat per-chrom H5s within each fold, align by TSV variant order, then average across folds."
    )
    p.add_argument("--fold-glob",           required=True)
    p.add_argument("--tsv-path",            required=True,
                   help="TSV used as input to ChromBPNet scoring, with chr/pos/variant_id columns")
    p.add_argument("--per-fold-file-glob",  default="**/*.h5")
    p.add_argument("--out",                 required=True)
    p.add_argument("--keep-keys",           default="")
    p.add_argument("--dataset-prefix",      default="")
    p.add_argument("--dry-run",             action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    keep_keys = None
    if args.keep_keys:
        keep_keys = set(s for s in args.keep_keys.split(",") if s)
    dataset_prefix = args.dataset_prefix if args.dataset_prefix else None

    print("[ARGS]", json.dumps({
        "fold_glob":          args.fold_glob,
        "tsv_path":           args.tsv_path,
        "per_fold_file_glob": args.per_fold_file_glob,
        "out":                args.out,
        "keep_keys":          sorted(list(keep_keys)) if keep_keys else None,
        "dataset_prefix":     dataset_prefix,
        "dry_run":            args.dry_run,
    }, indent=2))

    concat_then_average_folds(
        fold_glob=args.fold_glob,
        tsv_path=args.tsv_path,
        per_fold_file_glob=args.per_fold_file_glob,
        out_path=args.out,
        keep_keys=keep_keys,
        dataset_prefix=dataset_prefix,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())