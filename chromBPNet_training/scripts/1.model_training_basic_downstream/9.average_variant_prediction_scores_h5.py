#!/usr/bin/env python3
# concat_then_average.py

import os, sys, glob, json, argparse
import h5py
import numpy as np
import hdf5plugin  # registers filters if present


def _is_same_shape_except_axis(a_shape, b_shape, axis):
    if len(a_shape) != len(b_shape):
        return False
    return all((i == axis) or (a_shape[i] == b_shape[i]) for i in range(len(a_shape)))


def _concat_variable_axis(arrays, key):
    """
    Concatenate arrays along the single axis whose size varies across inputs.
    All other axes must match. Raises if 0 or >1 axes vary.
    """
    if not arrays:
        raise ValueError(f"No arrays to concatenate for key {key}")
    if len(arrays) == 1:
        return arrays[0]

    shapes = [a.shape for a in arrays]
    rank = len(shapes[0])
    if any(len(s) != rank for s in shapes):
        raise ValueError(f"Rank mismatch for {key}: {shapes}")

    vary_axes = []
    for ax in range(rank):
        sizes = {s[ax] for s in shapes}
        if len(sizes) > 1:
            vary_axes.append(ax)
        else:
            # double-check equality across arrays for fixed axes
            if any(s[ax] != shapes[0][ax] for s in shapes[1:]):
                raise ValueError(f"Inconsistent fixed axis {ax} for {key}: {shapes}")

    if len(vary_axes) != 1:
        raise ValueError(f"Expected exactly one varying axis for {key}, found {vary_axes} with shapes {shapes}")

    axis = vary_axes[0]
    return np.concatenate(arrays, axis=axis)


def _read_all_datasets(h5, keep_keys=None, dataset_prefix=None):
    """
    Return {'group/dataset': ndarray} for all leaf datasets.
    - If keep_keys is provided, only load those exact keys.
    - Else if dataset_prefix is provided (e.g., 'observed'), only load keys starting with 'observed/'.
    """
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
        # default: include everything
        out[name] = obj[()]

    h5.visititems(visitor)
    return out


def _group_by_key(list_of_dicts):
    out = {}
    for d in list_of_dicts:
        for k, v in d.items():
            out.setdefault(k, []).append(v)
    return out


def concat_chroms_within_fold(per_chrom_h5_paths, keys_1d=(), keep_keys=None, dataset_prefix=None):
    per_file = []
    for p in sorted(per_chrom_h5_paths):
        with h5py.File(p, 'r') as f:
            per_file.append(_read_all_datasets(f, keep_keys=keep_keys, dataset_prefix=dataset_prefix))

    by_key = _group_by_key(per_file)
    fold_cat = {}

    for key, arrs in by_key.items():
        if key in keys_1d:
            for a in arrs:
                if a.ndim != 1:
                    raise ValueError(f"{key} expected 1D but saw shape {a.shape}")
            fold_cat[key] = np.concatenate(arrs, axis=0)
        else:
            # OLD: fold_cat[key] = _concat_last_axis(arrs, key)
            fold_cat[key] = _concat_variable_axis(arrs, key)

    return fold_cat



def _reindex_last_axis(arr, idx):
    slicer = [slice(None)] * (arr.ndim - 1) + [idx]
    return arr[tuple(slicer)]


def align_folds_by_ids(concatenated_per_fold, ids_key):
    ref_ids = concatenated_per_fold[0].get(ids_key, None)
    if ref_ids is None:
        raise KeyError(f"ids_key '{ids_key}' not found in first fold")

    ref_ids = np.asarray(ref_ids).astype(str)
    ref_map = {v: i for i, v in enumerate(ref_ids)}

    aligned = []
    for fi, d in enumerate(concatenated_per_fold):
        ids = d.get(ids_key, None)
        if ids is None:
            raise KeyError(f"ids_key '{ids_key}' not found in fold {fi}")
        ids = np.asarray(ids).astype(str)

        try:
            idx = np.array([ref_map[v] for v in ids], dtype=int)
        except KeyError as e:
            raise ValueError(f"Fold {fi} contains id not present in reference: {e.args[0]}")

        if len(ref_ids) != len(ids) or set(ref_ids) != set(ids):
            raise ValueError(f"ID sets differ: ref={len(ref_ids)} vs fold{fi}={len(ids)}")

        d_aligned = {}
        for k, arr in d.items():
            if arr.ndim >= 1 and arr.shape[-1] == ids.shape[0]:
                d_aligned[k] = _reindex_last_axis(arr, np.argsort(idx))
            else:
                d_aligned[k] = arr
        aligned.append(d_aligned)

    return aligned


def average_folds(concatenated_per_fold, cast_back_float16=True):
    keys = set(concatenated_per_fold[0].keys())
    for d in concatenated_per_fold[1:]:
        if set(d.keys()) != keys:
            missing = keys.symmetric_difference(set(d.keys()))
            raise ValueError(f"Key mismatch across folds; differing keys: {missing}")

    averages = {}
    for key in sorted(keys):
        stacks = [d[key] for d in concatenated_per_fold]
        shape0 = stacks[0].shape
        if any(s.shape != shape0 for s in stacks[1:]):
            raise ValueError(f"Shape mismatch across folds for key {key}: "
                             f"{[s.shape for s in stacks]}")
        dtype0 = stacks[0].dtype
        acc = np.zeros(shape0, dtype=(np.float32 if dtype0 == np.float16 else dtype0))
        for s in stacks:
            acc += s.astype(acc.dtype, copy=False)
        avg = acc / float(len(stacks))
        if cast_back_float16 and dtype0 == np.float16:
            avg = avg.astype(np.float16)
        averages[key] = avg
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
    per_fold_file_glob="**/*.h5",
    out_path="averaged_folds.h5",
    keys_1d=(),
    keep_keys=None,
    dataset_prefix=None,
    align_by_ids=False,
    variant_ids_key='variant_ids',
    dry_run=False,
):
    fold_dirs = sorted(glob.glob(fold_glob))
    if not fold_dirs:
        raise FileNotFoundError(f"No fold dirs matched: {fold_glob}")

    per_fold_concat = []
    for fd in fold_dirs:
        h5s = sorted(glob.glob(os.path.join(fd, per_fold_file_glob), recursive=True))
        if not h5s:
            print(f"[WARN] No H5 files in {fd} (pattern {per_fold_file_glob})")
            continue
        print(f"[INFO] Fold {fd}: {len(h5s)} per-chrom files")
        if dry_run:
            for p in h5s:
                print("       ", p)
        fold_cat = concat_chroms_within_fold(
            h5s, keys_1d=tuple(keys_1d), keep_keys=keep_keys, dataset_prefix=dataset_prefix
        )
        per_fold_concat.append(fold_cat)

    if not per_fold_concat:
        raise RuntimeError("No per-fold data assembled; check your patterns.")

    if align_by_ids:
        print(f"[INFO] Aligning folds by ids key: {variant_ids_key}")
        per_fold_concat = align_folds_by_ids(per_fold_concat, ids_key=variant_ids_key)

    print("[INFO] Averaging across folds")
    averaged = average_folds(per_fold_concat, cast_back_float16=True)

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
        description="Concatenate per-chrom H5s within each fold, then average arrays across folds."
    )
    p.add_argument("--fold-glob", required=True,
                   help="Glob for fold directories, e.g. /.../fold_*")
    p.add_argument("--per-fold-file-glob", default="**/*.h5",
                   help="Glob inside each fold to find per-chrom files (recursive allowed)")
    p.add_argument("--out", required=True, help="Output H5 path")
    p.add_argument("--keys-1d", default="",
                   help="Comma-separated keys that are 1D and should concat axis=0 (e.g., variant_ids)")
    p.add_argument("--keep-keys", default="",
                   help="Comma-separated full dataset keys to keep "
                        "(e.g. 'observed/allele1_pred_profiles,observed/allele2_pred_profiles'). "
                        "Leave empty to include all datasets (or combine with --dataset-prefix).")
    p.add_argument("--dataset-prefix", default="",
                   help="Only include datasets whose key starts with this prefix (e.g., 'observed')")
    p.add_argument("--align-by-ids", action="store_true",
                   help="Align all folds by the variant_ids key before averaging")
    p.add_argument("--variant-ids-key", default="variant_ids",
                   help="Dataset key holding variant ids (used with --align-by-ids)")
    p.add_argument("--dry-run", action="store_true",
                   help="Parse and report shapes without writing output")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    keys_1d = tuple([s for s in args.keys_1d.split(",") if s]) if args.keys_1d else tuple()
    keep_keys = None
    if args.keep_keys:
        keep_keys = set(s for s in args.keep_keys.split(",") if s)
    dataset_prefix = args.dataset_prefix if args.dataset_prefix else None

    print("[ARGS]", json.dumps({
        "fold_glob": args.fold_glob,
        "per_fold_file_glob": args.per_fold_file_glob,
        "out": args.out,
        "keys_1d": keys_1d,
        "keep_keys": sorted(list(keep_keys)) if keep_keys else None,
        "dataset_prefix": dataset_prefix,
        "align_by_ids": args.align_by_ids,
        "variant_ids_key": args.variant_ids_key,
        "dry_run": args.dry_run
    }, indent=2))

    concat_then_average_folds(
        fold_glob=args.fold_glob,
        per_fold_file_glob=args.per_fold_file_glob,
        out_path=args.out,
        keys_1d=keys_1d,
        keep_keys=keep_keys,
        dataset_prefix=dataset_prefix,
        align_by_ids=args.align_by_ids,
        variant_ids_key=args.variant_ids_key,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())


