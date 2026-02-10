#!/usr/bin/env python3

import pandas as pd
import numpy as np
import yaml
import argparse
from pathlib import Path
from typing import Optional, Tuple, Dict, List


def load_chr(res_pre: str, fin_pre: str, chr_n: int, min_p2p3: Optional[float] = 0.5):
    """
    Load Findr results + finemapping results for one chromosome.
    NOTE: min_p2p3 is OPTIONAL; set to None to keep everything.
    """
    print(f"\n\tProcessing chromosome {chr_n}...")

    res = pd.read_csv(f"{res_pre}{chr_n}.tsv", sep="\t")
    if (min_p2p3 is not None) and ("p2p3" in res.columns):
        res = res.loc[res["p2p3"] > float(min_p2p3)].copy()
    res["chromosome"] = chr_n

    fin = pd.read_csv(f"{fin_pre}{chr_n}.tsv", sep="\t")
    fin["chromosome"] = chr_n

    return res, fin


def pp_cutoff_for_fdr(res: pd.DataFrame, thd: float) -> float:
    """
    Scan p2p3 thresholds s in [0.5, 1) and estimate:
      FDR(s) = 1 - mean(p2p3 | p2p3 >= s)
    then choose the s whose FDR is closest to thd.

    Returns the selected p2p3 cutoff. If can't compute, returns +inf.
    """
    if res.empty or "p2p3" not in res.columns:
        return float("inf")

    steps = np.arange(0.5, 1.0, 0.001)
    fdr_l = []
    p = res["p2p3"]

    for s in steps:
        m = p.loc[p >= s].mean()
        fdr_l.append(1.0 - m if pd.notnull(m) else np.nan)

    fdr_scan = pd.DataFrame({"p2p3": steps, "FDR": fdr_l}).dropna()
    if fdr_scan.empty:
        return float("inf")

    sel = fdr_scan.iloc[(fdr_scan["FDR"] - thd).abs().argsort()[:1]]["p2p3"].values[0]
    return float(sel)


def add_fdr_labels(res: pd.DataFrame, thds: List[float]) -> Tuple[pd.DataFrame, Dict[float, float]]:
    """
    Add boolean columns fdr_sig_{int(thd*100)} to res based on computed cutoffs.
    Also return dict mapping thd -> cutoff used.
    """
    cutoffs: Dict[float, float] = {}
    out = res.copy()

    for thd in thds:
        cutoff = pp_cutoff_for_fdr(out, thd)
        cutoffs[thd] = cutoff

        sig_col = f"fdr_sig_{int(thd * 100)}"
        out[sig_col] = (out["p2p3"] >= cutoff) if np.isfinite(cutoff) else False

        # Provenance: store the cutoff used
        out[f"p2p3_cutoff_fdr_{int(thd * 100)}"] = cutoff

    return out, cutoffs


def annotate_and_write_chr(
    atype: str,
    res: pd.DataFrame,
    fin: pd.DataFrame,
    tss_fn: str,
    dist: int,
    thds: List[float],
    out_pre: str,
    chr_n: int,
):
    """
    Annotate results (coloc_triplet, cis/trans proximity, FDR labels) and write ONE file per chromosome.
    Keeps CIS only (no FDR filtering; FDR is just labeled).
    """
    out_a = {
        "chromatin": "caqtl_CD4T_combined_sig_shared_inpeak_findr_results_mediation_raw_annotated_",
        "gene": "eqtl_all_CD4T_cells_gene_expression_mediation_raw_annotated_",
    }

    Path(out_pre).mkdir(parents=True, exist_ok=True)
    out_fn = f"{out_pre}{out_a[atype]}chr{chr_n}.tsv"

    if res.empty:
        print(f"\t  (chr{chr_n}) Findr results empty after initial filters; writing empty annotated file.")
        res.to_csv(out_fn, sep="\t", index=False)
        print(f"\t  Wrote: {out_fn}")
        return

    # ---- Triplet info from fine-mapping table ----
    fin_filt = fin[["peak", "variant_id", "gene", "finemapped_cs_coloc"]].copy()
    fin_filt.columns = ["peak_name", "exposure_variable", "gene_name", "finemapped_cs_coloc"]
    fin_filt["coloc_triplet"] = True

    merged = pd.merge(
        res,
        fin_filt,
        on=["peak_name", "exposure_variable", "gene_name", "finemapped_cs_coloc"],
        how="left",
    )
    merged["coloc_triplet"] = merged["coloc_triplet"].fillna(False).astype(bool)

    # ---- Load TSS/feature positions ----
    tss = pd.read_csv(tss_fn, sep="\t", usecols=[0, 1, 3])
    tss["#chr"] = tss["#chr"].astype(str).str.replace("chr", "", regex=False).astype("int64")

    # ---- Add FDR label columns (no filtering) ----
    merged, _ = add_fdr_labels(merged, thds)

    # ---- Add proximity labels (cis/trans), then FILTER CIS ONLY ----
    merged["variant_position"] = merged["exposure_variable"].str.extract(r":(\d+)\[")[0].astype(int)

    if atype == "chromatin":
        merged = pd.merge(merged, tss, on="gene_name", how="left")
        merged = merged.dropna(subset=["#chr", "start"])

        merged["gene_ev_proximity"] = (
            (merged["#chr"] == merged["chromosome"])
            & (merged["variant_position"].between(merged["start"] - dist, merged["start"] + dist))
        ).map({True: "cis", False: "trans"})

        merged = merged.loc[merged["gene_ev_proximity"] == "cis"].copy()

        base_cols = [
            "finemapped_cs_coloc",
            "exposure_variable",
            "chromosome",
            "peak_name",
            "gene_name",
            "coloc_triplet",
            "gene_ev_proximity",
            "p2",
            "p3",
            "p4",
            "p5",
            "p2p3",
        ]

    elif atype == "gene":
        merged = pd.merge(merged, tss, on="peak_name", how="left")
        merged = merged.dropna(subset=["#chr", "start"])

        merged["peak_ev_proximity"] = (
            (merged["#chr"] == merged["chromosome"])
            & (merged["variant_position"].between(merged["start"] - dist, merged["start"] + dist))
        ).map({True: "cis", False: "trans"})

        merged = merged.loc[merged["peak_ev_proximity"] == "cis"].copy()

        base_cols = [
            "finemapped_cs_coloc",
            "exposure_variable",
            "chromosome",
            "gene_name",
            "peak_name",
            "coloc_triplet",
            "peak_ev_proximity",
            "p2",
            "p3",
            "p4",
            "p5",
            "p2p3",
        ]
    else:
        raise ValueError(f"Invalid atype '{atype}'. Expected 'chromatin' or 'gene'.")

    # Keep all FDR columns
    fdr_cols = [c for c in merged.columns if c.startswith("fdr_sig_") or c.startswith("p2p3_cutoff_fdr_")]
    annotated = merged[base_cols + fdr_cols].copy()

    annotated.to_csv(out_fn, sep="\t", index=False)
    print(f"\t  Wrote: {out_fn}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--atype", required=True, type=str, help="A-variable type: 'gene' or 'chromatin'")
    parser.add_argument(
        "--config",
        default="/gchm/cd4_qtl_paper_figures/figure_3/analysis/mediation_dec25/med_res_filtering/mediation_results_filter.yml",
        help="Path to config yml",
    )
    parser.add_argument(
        "--min_p2p3",
        type=float,
        default=0.5,
        help="Optional pre-filter on p2p3 before annotation (set e.g. 0.0 to keep most; or use --keep_all).",
    )
    parser.add_argument(
        "--keep_all",
        action="store_true",
        help="Keep all rows (overrides --min_p2p3).",
    )
    args = parser.parse_args()

    atype = args.atype
    if atype not in ("gene", "chromatin"):
        raise ValueError(f"Unsupported atype: {atype}. Expected 'gene' or 'chromatin'")

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    cpaths = config["paths"]
    cparms = config["parms"]

    dparms_full = {
        "chromatin": ["cres_pre", "cfin_pre", "mrna_fn"],
        "gene": ["eres_pre", "efin_pre", "peak_fn"],
    }
    res_key, fin_key, tss_key = dparms_full[atype]

    thds = [float(p["thd"]) for p in cparms.get("test_parms", [])]
    if not thds:
        raise ValueError("No 'test_parms' thresholds found in config (parms:test_parms).")

    dist = int(cparms["dist_parms"]["dist"])

    print(f"\nA-variable: {atype}")
    print(f"Thresholds (label only, no filtering): {thds}")
    print(f"Distance for cis/trans: {dist}")

    min_p2p3: Optional[float] = None if args.keep_all else float(args.min_p2p3)

    print("\nAnnotating per chromosome...")
    for chr_n in range(1, 23):
        res_chr, fin_chr = load_chr(cpaths[res_key], cpaths[fin_key], chr_n, min_p2p3=min_p2p3)
        annotate_and_write_chr(
            atype=atype,
            res=res_chr,
            fin=fin_chr,
            tss_fn=cpaths[tss_key],
            dist=dist,
            thds=thds,
            out_pre=cpaths["out_pre"],
            chr_n=chr_n,
        )

    print("\nComplete.\n")


if __name__ == "__main__":
    main()
