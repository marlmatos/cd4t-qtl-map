# QTL fine-mapping and colocalization

SuSiE fine-mapping of caQTL and eQTL summary statistics (from [`qtl_mapping/`](../qtl_mapping)), followed by colocalization: **eQTL vs caQTL**, and **each QTL type vs immune-related GWAS**.

## Pipeline

Steps are numbered in run order; each numbered R script has a matching SLURM submission `.sh` wrapper.

### 00 — Orchestration wrappers (per QTL type)
Loop over chromosomes 1–22, calling steps 01 → 02 → 03 for each region:
- `00.Run_finemap_pipeline_caQTL.sh` — caQTL fine-mapping (ATAC-seq peaks)
- `00.Run_finemap_pipeline_eQTL.sh` — eQTL fine-mapping (all-cell pseudobulk)

### 01 — Define fine-mapping regions
- `01.Define_finemap_regions.caQTL.R` / `01.Define_finemap_regions.eQTL.R` — from the independent (conditionally significant) cis-QTL results, define a genomic window around each genome-wide-significant SNP to fine-map

### 02 — LD calculation
- `02.Generate_LD.sh` — PLINK LD (`--r2 square`) for each region, from cohort genotypes (duplicate variants excluded)

### 03 — SuSiE fine-mapping
- `03.submit_caQTL_finemap.sh` / `03.submit_finemap.sh` — SLURM wrappers
- `03.susie_finemap_caQTLs.R` / `03.susie_finemap_eQTLs.R` — SuSiE fine-mapping of tensorQTL summary statistics per region; outputs credible sets

### 04 — GWAS acquisition and liftover
- `04.Download_GWAS_sumstats.sh` — download GWAS summary statistics (GWAS Catalog harmonised files, immune-trait study list)
- `04.Run_liftover.sh` / `04.liftover_finemap_QTLs.sh` — liftover fine-mapped QTL credible sets hg38 → hg19 (bedtools, SLURM array job per chromosome)
- `04.liftover_finemapped_QTLs.R` — companion script converting/annotating the lifted coordinates against dbSNP

### 05 — LD (1000 Genomes) and QTL–GWAS colocalization
- `05.Generate_LD_1000genomes.sh` — PLINK2 LD from 1000 Genomes phase 3 (GRCh38, EUR) for GWAS regions
- `05.Run_coloc.R` — generic `coloc`/SuSiE colocalization runner (single GWAS vs QTL)
- `05.Run_coloc_immune_GWAS.R` — eQTL vs immune GWAS colocalization
- `05.Run_coloc_immune_caQTL_GWAS.R` — caQTL vs immune GWAS colocalization
- `05.Submit_coloc.sh` / `05.Submit_coloc_caQTL.sh` — submit one coloc job per GWAS study (reads `preprocessed_folder_names.txt`)

### 06 — eQTL–caQTL colocalization
- `06.Run_ca_eqtl_coloc.sh` — SLURM array job (chr 1–22)
- `06.coloc_eQTL_caQTL.R` — `coloc` between fine-mapped caQTL and eQTL signals

### 07–08 — Downstream summaries
- `07.make_locus_plots.R` — regional association/locus plots for significant QTL–GWAS coloc hits
- `08.credible_sets_size.R` — summarize fine-mapped region widths and credible-set sizes

## Supporting / utility scripts
- `Make_plots.R` — summary plots across eQTL–GWAS coloc results, grouped by study year/trait
- `compare_independent_snps.R` — compare tensorQTL's independent significant SNPs against the fine-mapped credible sets
- `process_coloc_results.R` — merge/annotate coloc result tables against GWAS study metadata
- `prop_GWAS_with_QTL.R` — proportion of GWAS credible sets with a colocalizing QTL
- `finemap_ca_only.sh` / `finemap_only.sh` — standalone loop wrappers submitting `03.submit_caQTL_finemap.sh` / `03.submit_finemap.sh` across chromosomes 1–22
- `convert_to_plink.sh` / `download_data.sh` / `liftover_vcf.sh` — one-off reference-data prep (PLINK conversion, GWAS download, 1000 Genomes VCF liftover)

## Paths

These scripts contain absolute paths to lab/cluster storage, masked here as `/gchm` (personal scratch), `/gchm_lab/collab` (collaborator's shared project directory), and `/gchm_lab/collab_gwas` (shared GWAS fine-mapping output). Update these to your own storage locations before running.
