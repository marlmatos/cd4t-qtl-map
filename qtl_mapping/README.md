# QTL mapping

Genotype preparation and QTL mapping (tensorQTL) for chromatin accessibility (caQTL) and gene expression (eQTL) in CD4+ T cells.

## Directory structure

```
qtl_mapping/
├── genotypes/       # VCF QC, allele splitting, PC filtering (shared by caQTL + eQTL)
├── caqtl_mapping/    # ATAC-seq peak caQTL mapping
└── eqtl_mapping/     # scRNA-seq pseudobulk eQTL mapping
```

## `genotypes/`
Genotype QC shared by both QTL workflows:
- `001_filter_qc.sh` — per-chromosome VCF filtering/QC
- `002_split_alleles.sh` — split multiallelic sites
- `003_variant_filtering_pcs_MAF5.ipynb` — MAF filtering and genotype PC calculation

## `caqtl_mapping/`
Numbered pipeline (run in order) mapping caQTLs from ATAC-seq peaks with tensorQTL:
1. `001_prepare_peak_file` — build the ATAC-seq peak quantification matrix
2. `002_Identify_feature_PCs` — PCA on the peak matrix
3. `003_Process_input_peaks` — format peak input for tensorQTL
4. `004_Prepare_Genotypes_PCs` — genotype PC preparation
5. `005_Change_var_names` — standardize variant naming
6. `006_Splitting_genotyped_by_chr` — split PLINK files by chromosome
7. `007_Prepare_CovsFile` — build the covariates file
8. `008_Run_caQTL_narrowpeaks_tensor` / `008_caQTL_narrowpeaks_tensor-1mb.py` — run cis-caQTL mapping (±1Mb window, MAF>5%)

Each `.R`/`.py` script has a matching `.sh` SLURM submission wrapper.

## `eqtl_mapping/`
Numbered pipeline (run in order) mapping cis-eQTLs from scRNA-seq pseudobulk with tensorQTL:
1. `001_preparing_seurat_obj_for_QTL` — subset Seurat object to CD4+ cells
2. `002_GetPseudobulk_whole` — compute per-sample pseudobulk expression
3. `003_IdentifyPCs` — PCA on expression (requires `NormalizePseudobulk.R` in this directory — not yet committed, see below)
4. `004_ProcessExpression` — format expression input for tensorQTL
5. `005_PreparecovsFile` / `005_PreparecovsFile_ncells` — build covariates file(s)
6. `006_ReFiltering_genotypes` — re-filter genotypes for this cohort
7. `007_split_plinkfiles_by_chr` — split PLINK files by chromosome
8. `008_ciseQTLTensor_allcells` — run cis-eQTL mapping (±1Mb window, MAF>5%, read-depth covariate)

`plotting_pcs.ipynb` and `003_common_samples_gfgex_wgs.csv` are supporting QC/metadata files.

**Note:** `003_IdentifyPCs.R` sources `NormalizePseudobulk.R`, a pseudobulk-normalization utility that has not yet been added to this directory — add it at `qtl_mapping/eqtl_mapping/NormalizePseudobulk.R` before running this step.
