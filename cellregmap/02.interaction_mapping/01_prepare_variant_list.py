import os
import pandas as pd
import csv

# -----------------------------
# USER INPUTS
# -----------------------------

res = "34"

finemap_dir = "gcgl/sghatan/marlis_pj/coloc/SuSiE_finemap_credible_sets/All_CD4T_cells"

# Example pattern:
# All_CD4T_cells_chr1_credible_sets.txt
finemap_pattern = "All_CD4T_cells_chr{chrom}_credible_sets.txt"

nominal_dir = "~/cd4_QTL_analysis/03_Run_cisQTL_perchr/analysis/cis_eQTLs/all_CD4T_cells_MAF5"

# Replace this with your actual file pattern
# Example possibilities:
# "All_CD4T_cells_chr{chrom}.cis_qtl_pairs.txt.gz"
# "allcells.cis_qtl_pairs.chr{chrom}.csv"
nominal_pattern = "allcells.cis_qtl_pairs.chr{chrom}.csv"

outdir = "~/cd4_CellRegMap/002_interaction_analysis/results/results_042426"
os.makedirs(outdir, exist_ok=True)


# -----------------------------
# READ FINEMAPPING CREDIBLE SET FILES
# -----------------------------

finemap_list = []

for chrom in range(1, 23):
    finemap_path = os.path.join(
        finemap_dir,
        finemap_pattern.format(chrom=chrom)
    )

    if not os.path.exists(finemap_path):
        print(f"Warning: fine-mapping file not found for chr{chrom}: {finemap_path}")
        continue

    try:
        df = pd.read_csv(finemap_path, sep="\t")

        if df.empty:
            print(f"Warning: fine-mapping file for chr{chrom} is empty.")
            continue

        df["chrom_file"] = chrom
        finemap_list.append(df)

    except pd.errors.EmptyDataError:
        print(f"Warning: fine-mapping file for chr{chrom} is empty. Skipping.")


if not finemap_list:
    raise FileNotFoundError("No fine-mapping credible set files were found.")

finemap_df = pd.concat(finemap_list, ignore_index=True)

print("Fine-mapping dataframe shape:", finemap_df.shape)
print(finemap_df.columns.tolist())


# -----------------------------
# CLEAN AND DEFINE CREDIBLE SETS
# -----------------------------

required_finemap_cols = ["variant_id", "gene", "region", "cs", "pip"]

missing_cols = [col for col in required_finemap_cols if col not in finemap_df.columns]
if missing_cols:
    raise ValueError(f"Missing required fine-mapping columns: {missing_cols}")

finemap_df = finemap_df.copy()

finemap_df["variant_id"] = finemap_df["variant_id"].astype(str).str.strip()
finemap_df["gene"] = finemap_df["gene"].astype(str).str.strip()
finemap_df["region"] = finemap_df["region"].astype(str).str.strip()
finemap_df["cs"] = finemap_df["cs"].astype(str).str.strip()
finemap_df["pip"] = pd.to_numeric(finemap_df["pip"], errors="coerce")

finemap_df = finemap_df.dropna(
    subset=["variant_id", "gene", "region", "cs", "pip"]
)

# Credible set = region + gene + cs
finemap_df["credible_set_id"] = (
    finemap_df["region"] + "__" +
    finemap_df["gene"] + "__CS" +
    finemap_df["cs"]
)

# -----------------------------
# GET LEAD VARIANT PER CREDIBLE SET
# -----------------------------

lead_idx = (
    finemap_df
    .groupby(["region", "gene", "cs"])["pip"]
    .idxmax()
)

lead_finemap = (
    finemap_df
    .loc[lead_idx]
    .copy()
    .sort_values(["gene", "region", "cs"])
)

lead_finemap = lead_finemap.rename(
    columns={
        "gene": "phenotype_id",
        "variant_id": "lead_variant_id",
        "pip": "lead_pip"
    }
)

print("Number of fine-mapped credible sets:", lead_finemap.shape[0])
print("Number of unique lead variants:", lead_finemap["lead_variant_id"].nunique())
print("Number of unique genes:", lead_finemap["phenotype_id"].nunique())

# -----------------------------
# READ NOMINAL CIS-eQTL ASSOCIATION FILES
# -----------------------------

nominal_list = []

for chrom in range(1, 23):
    nominal_path = os.path.join(
        nominal_dir,
        nominal_pattern.format(chrom=chrom)
    )

    if not os.path.exists(nominal_path):
        print(f"Warning: nominal association file not found for chr{chrom}: {nominal_path}")
        continue

    try:
        # Change sep if your nominal files are comma-separated
        df = pd.read_csv(nominal_path, sep="\t")

        if df.empty:
            print(f"Warning: nominal association file for chr{chrom} is empty.")
            continue

        df["chrom"] = chrom
        nominal_list.append(df)

    except pd.errors.EmptyDataError:
        print(f"Warning: nominal association file for chr{chrom} is empty. Skipping.")


if not nominal_list:
    raise FileNotFoundError("No nominal cis-eQTL association files were found.")

nominal_df = pd.concat(nominal_list, ignore_index=True)

print("Nominal association dataframe shape:", nominal_df.shape)
print(nominal_df.columns.tolist())

# -----------------------------
# CLEAN NOMINAL ASSOCIATION DATA
# -----------------------------

if "phenotype_id" not in nominal_df.columns:
    if "gene" in nominal_df.columns:
        nominal_df = nominal_df.rename(columns={"gene": "phenotype_id"})
    else:
        raise ValueError("Nominal association files need either 'phenotype_id' or 'gene' column.")

required_nominal_cols = ["variant_id", "phenotype_id", "pval_nominal"]

missing_cols = [col for col in required_nominal_cols if col not in nominal_df.columns]
if missing_cols:
    raise ValueError(f"Missing required nominal association columns: {missing_cols}")

nominal_df = nominal_df.copy()

nominal_df["variant_id"] = nominal_df["variant_id"].astype(str).str.strip()
nominal_df["phenotype_id"] = nominal_df["phenotype_id"].astype(str).str.strip()
nominal_df["pval_nominal"] = pd.to_numeric(nominal_df["pval_nominal"], errors="coerce")

# Optional: keep Cell.type if present
nominal_keep_cols = ["variant_id", "phenotype_id", "chrom", "pval_nominal"]

if "Cell.type" in nominal_df.columns:
    nominal_keep_cols.append("Cell.type")

# Keep only the columns needed for CellRegMap input
nominal_df = nominal_df[nominal_keep_cols].drop_duplicates()

# -----------------------------
# MERGE LEAD FINEMAPPED VARIANTS WITH NOMINAL STATS
# -----------------------------

cd4_eqtl = lead_finemap.merge(
    nominal_df,
    left_on=["lead_variant_id", "phenotype_id"],
    right_on=["variant_id", "phenotype_id"],
    how="left"
)

# Check merge success
n_missing_p = cd4_eqtl["pval_nominal"].isna().sum()

print("Merged CellRegMap input shape:", cd4_eqtl.shape)
print("Number of lead variants missing nominal p-value:", n_missing_p)

if n_missing_p > 0:
    print("Example missing nominal p-values:")
    print(
        cd4_eqtl.loc[
            cd4_eqtl["pval_nominal"].isna(),
            ["lead_variant_id", "phenotype_id", "region", "cs", "lead_pip"]
        ].head()
    )


# -----------------------------
# FINAL CLEANING BEFORE SELECTION
# -----------------------------

# Rename lead variant back to variant_id for compatibility with old downstream code
cd4_eqtl = cd4_eqtl.drop(columns=["variant_id"], errors="ignore")
cd4_eqtl = cd4_eqtl.rename(columns={"lead_variant_id": "variant_id"})

# Keep only rows with nominal p-values
cd4_eqtl = cd4_eqtl.dropna(subset=["pval_nominal"]).copy()

# Make sure required columns exist
if "Cell.type" not in cd4_eqtl.columns:
    cd4_eqtl["Cell.type"] = "allcells"

# Prefer nominal chrom column if it exists; otherwise derive from variant_id
if "chrom" not in cd4_eqtl.columns:
    cd4_eqtl["chrom"] = (
        cd4_eqtl["variant_id"]
        .str.split(":")
        .str[0]
        .str.replace("chr", "", regex=False)
        .astype(int)
    )

print("Merged fine-mapped lead table shape:", cd4_eqtl.shape)
print("Unique genes before one-SNP-per-gene filtering:", cd4_eqtl["phenotype_id"].nunique())
print("Unique lead variants before filtering:", cd4_eqtl["variant_id"].nunique())
print("Unique credible sets before filtering:", cd4_eqtl["credible_set_id"].nunique())

# -----------------------------
# OPTION A: ONE LEAD SNP PER GENE
# -----------------------------

# Keep the full version temporarily, with lead_pip and credible_set_id still present
cd4_eqtl_one_snp_per_gene_full = (
    cd4_eqtl
    .sort_values(
        ["phenotype_id", "lead_pip", "pval_nominal"],
        ascending=[True, False, True]
    )
    .drop_duplicates(subset=["phenotype_id"], keep="first")
    .copy()
)

print("Option A: one SNP per gene")
print("Number of genes:", cd4_eqtl_one_snp_per_gene_full["phenotype_id"].nunique())
print("Number of variants:", cd4_eqtl_one_snp_per_gene_full["variant_id"].nunique())
print("Number of rows:", cd4_eqtl_one_snp_per_gene_full.shape[0])

# -----------------------------
# FINAL CELLREGMAP INPUT FORMAT
# -----------------------------

final_cols = ["variant_id", "phenotype_id", "Cell.type", "chrom", "pval_nominal"]

cd4_eqtl_one_snp_per_gene = (
    cd4_eqtl_one_snp_per_gene_full[final_cols]
    .copy()
    .sort_values(["phenotype_id", "pval_nominal"])
)

# Enforce column dtypes/format
cd4_eqtl_one_snp_per_gene["Cell.type"] = (
    cd4_eqtl_one_snp_per_gene["Cell.type"]
    .fillna("allcells")
)

cd4_eqtl_one_snp_per_gene["chrom"] = (
    cd4_eqtl_one_snp_per_gene["chrom"]
    .astype(int)
)

cd4_eqtl_one_snp_per_gene["pval_nominal"] = pd.to_numeric(
    cd4_eqtl_one_snp_per_gene["pval_nominal"],
    errors="coerce"
)

print(cd4_eqtl_one_snp_per_gene.head())
print(cd4_eqtl_one_snp_per_gene.columns.tolist())

# -----------------------------
# SAVE OUTPUTS
# -----------------------------
outdir2="~/cd4_CellRegMap/002_interaction_analysis/results/results_042426"


cd4_eqtl_one_snp_per_gene.to_csv(
    os.path.join(outdir2, f"cis_eQTLS_allcells_one_lead_snp_per_gene_res{res}.csv"),
    index=False,
    quoting=csv.QUOTE_ALL
)

gene_list_one = (
    cd4_eqtl_one_snp_per_gene[["phenotype_id"]]
    .drop_duplicates()
    .sort_values("phenotype_id")
)

variant_list_one = (
    cd4_eqtl_one_snp_per_gene[["variant_id"]]
    .drop_duplicates()
    .sort_values("variant_id")
)

gene_list_one.to_csv(
    os.path.join(outdir2, f"egene_list_one_lead_snp_per_gene_res{res}.txt"),
    index=False,
    header=False,
    quoting=csv.QUOTE_ALL
)

variant_list_one.to_csv(
    os.path.join(outdir2, f"variant_list_one_lead_snp_per_gene_res{res}.txt"),
    index=False,
    header=False,
    quoting=csv.QUOTE_ALL
)


# -----------------------------
# OPTION B: MULTIPLE INDEPENDENT SNPs PER GENE
# -----------------------------

cd4_eqtl_independent = cd4_eqtl.copy()

# Count number of independent tests per gene
n_tests_per_gene = (
    cd4_eqtl_independent
    .groupby("phenotype_id")["credible_set_id"]
    .nunique()
    .reset_index(name="n_independent_tests_for_gene")
)

cd4_eqtl_independent = cd4_eqtl_independent.merge(
    n_tests_per_gene,
    on="phenotype_id",
    how="left"
)

# Bonferroni correction within gene
cd4_eqtl_independent["pval_nominal_bonf_gene"] = (
    cd4_eqtl_independent["pval_nominal"] *
    cd4_eqtl_independent["n_independent_tests_for_gene"]
).clip(upper=1.0)

# Sort for readability
cd4_eqtl_independent = cd4_eqtl_independent.sort_values(
    ["phenotype_id", "pval_nominal_bonf_gene", "pval_nominal"]
)

gene_list_independent = (
    cd4_eqtl_independent[["phenotype_id"]]
    .drop_duplicates()
    .sort_values("phenotype_id")
)

variant_list_independent = (
    cd4_eqtl_independent[["variant_id"]]
    .drop_duplicates()
    .sort_values("variant_id")
)

print("Option B: multiple independent SNPs per gene")
print("Number of genes:", cd4_eqtl_independent["phenotype_id"].nunique())
print("Number of variants:", cd4_eqtl_independent["variant_id"].nunique())
print("Number of independent SNP-gene tests:", cd4_eqtl_independent.shape[0])
print(
    "Genes with >1 independent SNP:",
    (cd4_eqtl_independent.groupby("phenotype_id").size() > 1).sum()
)

# -----------------------------
# FINAL CELLREGMAP INPUT FORMAT FOR OPTION B
# -----------------------------

final_cols = ["variant_id", "phenotype_id", "Cell.type", "chrom", "pval_nominal"]

cd4_eqtl_independent_cellregmap = (
    cd4_eqtl_independent[final_cols]
    .copy()
    .sort_values(["phenotype_id", "pval_nominal"])
)

# Enforce column dtypes/format
cd4_eqtl_independent_cellregmap["Cell.type"] = (
    cd4_eqtl_independent_cellregmap["Cell.type"]
    .fillna("allcells")
)

cd4_eqtl_independent_cellregmap["chrom"] = (
    cd4_eqtl_independent_cellregmap["chrom"]
    .astype(int)
)

cd4_eqtl_independent_cellregmap["pval_nominal"] = pd.to_numeric(
    cd4_eqtl_independent_cellregmap["pval_nominal"],
    errors="coerce"
)

print("Final Option B CellRegMap input:")
print(cd4_eqtl_independent_cellregmap.head())
print(cd4_eqtl_independent_cellregmap.columns.tolist())
print(cd4_eqtl_independent_cellregmap.shape)

# -----------------------------
# SAVE OPTION B OUTPUTS
# -----------------------------

# 1. CellRegMap-ready file with exact required columns only
cd4_eqtl_independent_cellregmap.to_csv(
    os.path.join(outdir2, f"cis_eQTLS_allcells_independent_lead_snps_per_gene_res{res}.csv"),
    index=False,
    quoting=csv.QUOTE_ALL
)

# 2. Full file with metadata, lead_pip, credible_set_id, and Bonferroni correction
cd4_eqtl_independent.to_csv(
    os.path.join(outdir2, f"cis_eQTLS_allcells_independent_lead_snps_per_gene_full_res{res}.csv"),
    index=False,
    quoting=csv.QUOTE_ALL
)

# 3. Gene list
gene_list_independent.to_csv(
    os.path.join(outdir2, f"egene_list_independent_lead_snps_per_gene_res{res}.txt"),
    index=False,
    header=False,
    quoting=csv.QUOTE_ALL
)

# 4. Variant list
variant_list_independent.to_csv(
    os.path.join(outdir2, f"variant_list_independent_lead_snps_per_gene_res{res}.txt"),
    index=False,
    header=False,
    quoting=csv.QUOTE_ALL
)

print("Saved Option B outputs:")
print(f"  CellRegMap input: cis_eQTLS_allcells_independent_lead_snps_per_gene_res{res}.csv")
print(f"  Full annotated file: cis_eQTLS_allcells_independent_lead_snps_per_gene_full_res{res}.csv")
print(f"  Gene list: egene_list_independent_lead_snps_per_gene_res{res}.txt")
print(f"  Variant list: variant_list_independent_lead_snps_per_gene_res{res}.txt")