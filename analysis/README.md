# Downstream integrative analysis

Integrative analyses combining caQTL, eQTL, fine-mapping, colocalization, chromBPNet predictions, and GWAS loci, used to generate manuscript results and figures.

## Directory structure

```
analysis/
├── 01-13...   # Numbered analysis scripts (run largely independently, some build on prior outputs)
├── helper_functions.R
├── utils/      # Shared plotting/annotation helper functions
└── resources/  # Cached reference annotations (e.g. ENCODE chromatin states)
```

## Local path configs

Scripts `01`, `02`, `03`, `11`, `12`, and `13` each source `paths_config.R` for cluster-specific absolute paths. That file is gitignored (`.gitignore`: `*/paths_config.R`) since it points at private data locations — copy [`paths_config.R.example`](paths_config.R.example) to `paths_config.R` in this directory and fill in your paths. Each script can then be run on its own, in any order. The same pattern applies to `plotting_notebooks/config.local.R` (see [`plotting_notebooks/README.md`](../plotting_notebooks/README.md)) and `chromBPNet_training/scripts/1.model_training_basic_downstream/utils/paths_config.R`.

## Scripts
1. `01.run_peak_coaccessibility_analysis.R` — co-accessibility between ATAC-seq peaks
2. `02.merging_finemapping_caqtl_eqtl_coloc_res.R` — integrate fine-mapped caQTL/eQTL credible sets with eQTL–caQTL colocalization results
3. `03.chromatin_tss_informed_epi_annotations.R` — annotate accessibility peaks by genomic context (promoter/distal) and ChromBPNet motif content
4. `04.findr_inference_results_filter.py` — filter Findr causal inference results
5. `05.findr_inference_cocat_filtered_results.sh` — concatenate filtered Findr results across chromosomes/batches
6. `06.findr_inference_locus_glm_enrichments.R` — Findr causal inference for mediated QTLs; locus-level GLM enrichment tests
7. `07.chrombpnet_qtl_variant_peakcoverage.R` — quantify accessibility signal in ChromBPNet-supported peaks
8. `08.TFmotif_hit_overlap_hitcaller.R` — overlap TF-MoDISco motif instances with hitcaller output
9. `09.TFmotif_hit_overlap_hitcaller_plus_neighbors.R` — same, including neighboring motif hits
10. `10.motif_hits_TF_family_enrich_promoter.R` — Fisher's exact test for TF-family motif enrichment in promoters vs. distal elements
11. `11.chrombpnet_var_liftover_hg38_to_hg19.R` — liftover ChromBPNet-scored variants from hg38 to hg19
12. `12.creating_summary_gwas_chrombpnet_coloc.R` — ChromBPNet support at GWAS and QTL–GWAS colocalized loci; builds Supplementary Table 8
13. `13.cbpnet_enrich_gwas_coloc_v2.R` — ChromBPNet prioritization enrichment among molQTL-colocalized GWAS credible-set variants

`helper_functions.R` — shared helpers (e.g. fine-mapping dataframe prep) sourced by the scripts above.

## `utils/`
- `metadata_helpers.R` — sample/cohort metadata utilities
- `color_pallete_helper.R` — shared plotting color palettes
- `track_plots_helpers.R` — genome track plotting helpers
- `qtl_violion_helper.R` — QTL violin plot helpers
- `dynamic_mofa_plot_beta_helper.R` — MOFA factor/beta plotting helpers
- `extract_variant_shap_helper.R` / `extract_variant_prediction_shap_helper.R` / `extract_variant_prediction_shap_helper_v2.R` — extract per-variant ChromBPNet SHAP/prediction scores
- `cbpnet_extract_allelic_contrib_plus_atrr_bigwigs.R` / `_v2.R` — extract allelic contribution scores and attribution bigWigs from ChromBPNet

## `resources/`
Cached reference annotations used by the scripts above, including ENCODE chromatin-state pages (`ENOCDE_chromatin_state/`) and an associated file index (`files.txt`).
