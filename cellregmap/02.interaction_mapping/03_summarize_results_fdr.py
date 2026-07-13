##########
## Script to summarize CellRegMap Interaction Results
## Author: Ana Cuomo
## Modified by: Marliette Matos
## Date: 08/05/2024
#########


#!/usr/bin/env python

import os
import glob
import argparse
import pandas as pd
import numpy as np
import re
from statsmodels.stats.multitest import fdrcorrection


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize CellRegMap interaction results and apply multiple testing correction."
    )

    parser.add_argument(
        "--res",
        type=int,
        required=True,
        help="Resolution value, e.g. 4."
    )

    parser.add_argument(
        "--results_path",
        type=str,
        required=True,
        help="Path to directory containing per-gene CellRegMap result CSV files."
    )

    parser.add_argument(
        "--out_results",
        type=str,
        required=True,
        help="Output directory where summary file will be written."
    )

    parser.add_argument(
        "--outfile",
        type=str,
        default=None,
        help="Output filename. Default: cellregmap_interaction_summary_res_<res>.csv"
    )

    parser.add_argument(
        "--correction",
        type=str,
        choices=["fdr", "bonferroni_fdr"],
        default="fdr",
        help=(
            "Multiple testing correction strategy. "
            "'fdr' = one-step FDR across all tests. "
            "'bonferroni_fdr' = Bonferroni within gene, then FDR."
        )
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="FDR significance threshold. Default: 0.05."
    )

    parser.add_argument(
        "--file_pattern",
        type=str,
        default="*.csv",
        help="File pattern for result files. Default: *.csv."
    )

    return parser.parse_args()


def get_gene_from_filename(file_path):
    basename = os.path.basename(file_path)

    pattern = r"^(.+?)_([^_]+)_([0-9]+)_b38_([^_]+)_([^_]+)_results\.csv$"
    match = re.match(pattern, basename)

    if match is None:
        raise ValueError(f"Could not parse gene from filename: {basename}")

    return match.group(1)


def read_cellregmap_file(file_path):
    """
    Read one CellRegMap result file and return a cleaned dataframe.
    """
    df = pd.read_csv(file_path, index_col=0)

    if df.shape[0] == 0:
        return None

    required_cols = {"pv", "variant", "chrom"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(
            f"File {file_path} is missing required columns: {missing}"
        )

    gene = get_gene_from_filename(file_path)

    out = pd.DataFrame({
        "gene": gene,
        "chrom": df["chrom"].values,
        "snpID": df["variant"].values,
        "pv_raw": pd.to_numeric(df["pv"], errors="coerce")
    })

    # Replace NA p-values with 1
    out["pv_raw"] = out["pv_raw"].fillna(1.0)

    # Keep p-values bounded
    out["pv_raw"] = out["pv_raw"].clip(lower=0.0, upper=1.0)

    return out


def apply_correction(df, correction="fdr", alpha=0.05):
    """
    Apply either:
    1. FDR only across all independent signal-gene tests
    2. Bonferroni within gene followed by FDR
    """

    df = df.copy()

    # Number of independent lead variants / signals tested per gene
    df["n_independent_signals_for_gene"] = (
        df.groupby("gene")["snpID"]
          .transform("nunique")
    )

    if correction == "fdr":
        df["pv_for_fdr"] = df["pv_raw"]

    elif correction == "bonferroni_fdr":
        df["pv_bonf_gene"] = (
            df["pv_raw"] * df["n_independent_signals_for_gene"]
        ).clip(upper=1.0)

        df["pv_for_fdr"] = df["pv_bonf_gene"]

    else:
        raise ValueError(f"Unknown correction mode: {correction}")

    reject, pv_adj = fdrcorrection(
        df["pv_for_fdr"].values,
        alpha=alpha,
        method="indep",
        is_sorted=False
    )

    df["pv_adj"] = pv_adj
    df["reject"] = reject
    df["correction"] = correction
    df["alpha"] = alpha

    df = df.sort_values("pv_adj")

    return df


def main():
    args = parse_args()

    os.makedirs(args.out_results, exist_ok=True)

    if args.outfile is None:
        args.outfile = f"cellregmap_interaction_summary_res_{args.res}_{args.correction}.csv"

    search_path = os.path.join(args.results_path, args.file_pattern)
    files = sorted(glob.glob(search_path))

    if len(files) == 0:
        raise FileNotFoundError(f"No files found matching: {search_path}")

    print(f"Found {len(files)} files.")
    print(f"Correction mode: {args.correction}")

    all_results = []
    count_success = 0
    count_failed = 0
    count_empty = 0

    for i, file_path in enumerate(files, start=1):
        if i % 500 == 0:
            print(f"Processed {i} files...")

        try:
            df_file = read_cellregmap_file(file_path)

            if df_file is None:
                count_empty += 1
                continue

            all_results.append(df_file)
            count_success += 1

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            count_failed += 1
            continue

    if len(all_results) == 0:
        raise RuntimeError("No valid result files were processed.")

    df = pd.concat(all_results, axis=0, ignore_index=True)

    print(f"Successfully processed files: {count_success}")
    print(f"Empty files skipped: {count_empty}")
    print(f"Failed files: {count_failed}")
    print(f"Total tests before correction: {df.shape[0]}")
    print(f"Total genes tested: {df['gene'].nunique()}")

    df = apply_correction(
        df,
        correction=args.correction,
        alpha=args.alpha
    )

    out_path = os.path.join(args.out_results, args.outfile)
    df.to_csv(out_path, index=False)

    print(f"Saved results to: {out_path}")
    print(f"Significant tests at FDR < {args.alpha}: {df['reject'].sum()}")


if __name__ == "__main__":
    main()