# Main manuscript figure notebooks

R Notebooks (`.Rmd`) that generate the main figures of the manuscript, along with their knitted HTML output.

## Contents
- `plotting_notebook_figure_1.Rmd` … `plotting_notebook_figure_6.Rmd` — one notebook per main figure (Figures 1–6)
- `plotting_notebook_figure_*.nb.html` — corresponding knitted HTML notebooks
- `knit_all.R` — utility to knit all `plotting_notebook_figure_*.Rmd` files in this directory (or a directory passed as an argument) in one run

Each notebook starts with `source("config.local.R")`, a gitignored file holding your local absolute paths (kept out of git for the same reason as `paths_config.R` — see [`analysis/README.md`](../analysis/README.md#local-path-configs)). Copy [`config.local.R.example`](config.local.R.example) to `config.local.R` in this directory and fill in your paths before knitting. Notebooks also source shared color-palette/plotting helpers from [`analysis/utils/`](../analysis/utils).
