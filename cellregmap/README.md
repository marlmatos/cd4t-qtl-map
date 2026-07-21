# CellRegMap interaction mapping

Single-cell gene-by-context (GxC) interaction eQTL mapping in CD4+ T cells using [CellRegMap](https://github.com/annacuomo/CellRegMap_analyses).

## Directory structure

```
cellregmap/
├── 01.preprocessing/        # Metacells, pseudobulk, MOFA phenotype/context preparation
└── 02.interaction_mapping/   # CellRegMap interaction testing, FDR summary, effect sizes
```

## `01.preprocessing/`
1. `001_convert_sobj_to_sce.R` — convert Seurat object to SCE / h5ad
2. `002_gex_cell_QC.py` — cell-level QC
3. `003_create_metacells.py` — build metacells (adapted from CellRegMap_analyses)
4. `004_pseudobulk_from_pseudocells.py` — pseudobulk expression from metacells
5. `005_prepare_input_file.py` — build the eQTL variant–gene pair test file
6. `006_prepare_phenotypes_MOFA.py` — prepare pseudobulk metacell input for MOFA
7. `007_extract_mofa_results.py` — extract MOFA factors used as GxC contexts

## `02.interaction_mapping/`
1. `01_prepare_variant_list` — build the variant list to test
2. `02_interaction_test*` — run CellRegMap interaction tests (genotype × MOFA-factor context), expanding genotype/kinship matrices across pseudocells; `_one_var_pergene` variant restricts to one variant per gene
3. `03_summarize_results_fdr` — aggregate results and apply multiple-testing correction (adapted from Ana Cuomo's CellRegMap_analyses)
4. `04_multiple_testing_qvalue.ipynb` — q-value diagnostics
5. `05_estimate_betas_sig_interactions` — estimate effect sizes for significant interactions

Each numbered Python script has a matching `.sh` SLURM submission wrapper.
