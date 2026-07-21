# chromBPNet training and variant scoring

Training, interpretation, and variant scoring of [chromBPNet](https://github.com/kundajelab/chrombpnet) models on CD4+ T cell ATAC-seq data.

## Directory structure

```
chromBPNet_training/
└── scripts/
    ├── 1.model_training_basic_downstream/  # Peak calling → fold splitting → bias/chromBPNet training → contribution/prediction scores
    │   └── utils/                          # Path config, variant annotation, cross-fold score summary
    ├── 2.denovo_TF_discovery/              # De novo motif discovery (TF-MoDISco) + reference motif sets
    │   └── references/
    ├── 3.variant_scoring_TF_hits/          # Variant effect prediction + SHAP contribution scores across folds
    ├── 4.TF_disruption_scores/             # TF motif hit calling on variants (finemo hitcaller)
    └── utils/                              # Locus/track plotting helpers
```

`utils/paths_config.R` (gitignored, cluster-specific) is required by the numbered pipeline below — copy [`utils/paths_config.R.example`](scripts/1.model_training_basic_downstream/utils/paths_config.R.example) to `utils/paths_config.R` and fill in your paths first.

## `1.model_training_basic_downstream/`
Numbered end-to-end pipeline (run in order), each with a SLURM `.sh` wrapper:
1. `1.select_samples.ipynb` — select high-read-depth ATAC-seq samples for training
2. `2.peak_calling.sh` — MACS3 peak calling
3. `3.split_folds_generate_nonpeaks.sh` — generate cross-validation folds and non-peak background regions
4. `4.train_bias_model_singularity.sh` — train the Tn5-bias model
5. `5.train_chrombpnet_model_singularity.sh` — train chromBPNet models (per fold)
6. `7.contribution_scores_singularity.sh` — compute per-base contribution (SHAP) scores
7. `8.prediction_scores_singularity.sh` — compute predicted accessibility signal
8. `9.average_contribution_scores_h5.py` / `9.run_average_contribution_scores_h5.sh` — average contribution scores across folds
9. `10.average_contribution_prediction_scores_bw.sh` — convert averaged scores to bigWig
10. `11.denovo_motif_discovery_singularity*.sh` — run TF-MoDISco on averaged contribution scores

## `2.denovo_TF_discovery/`
De novo motif discovery (TF-MoDISco) run with an expanded TF reference set (`references/`).

## `3.variant_scoring_TF_hits/`
Score candidate variants with trained chromBPNet models:
1. `1.prepare_variant_list.ipynb` — build the variant list to score
2. `2.variant_prediction_scores.sh` — predicted allelic effect scores
3. `3.average_across_folds.sh` / `4.average_variant_prediction_scores_h5.py` — average prediction scores across folds
4. `5.exploring_variant_contribution.ipynb` — exploratory QC of variant contribution scores
5. `6.calculate_variant_contribution_scores_shap.sh` — per-variant SHAP contribution scores
6. `7.average_variant_contribution_scores_h5.py` — average contribution scores across folds
7. `8.visualize_variant_shap.ipynb` — visualize variant-level SHAP tracks

## `4.TF_disruption_scores/`
`1.finemo_hicaller_FILTVARS.sh` — call TF motif hits disrupted by scored variants (finemo hitcaller).

## `utils/`
`locus_plotting_helper.R` / `locus_plotting_helper_extended.R` — shared helpers for locus/track plots used across the analyses above.
