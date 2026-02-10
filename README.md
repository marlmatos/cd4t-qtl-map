# CD4+ T Cell Multi-Omic QTL Mapping

Analysis code and figure-generation scripts for **"Integrating multi-omic QTLs and predictive models reveals regulatory architectures at immune-related GWAS loci in CD4+ T cells."**

Matos MR, et al. Integrating multi-omic QTLs and predictive models reveals regulatory architectures at immune related GWAS loci in CD4+ T cells. medRxiv. 2026;2026.01.27.26344979. doi:10.64898/2026.01.27.26344979.

## Overview

This repository contains comprehensive analysis pipelines for multi-omic QTL mapping in CD4+ T cells, integrating:
- **ATAC-seq** (chromatin accessibility)
- **scRNA-seq** (single-cell RNA sequencing)
- **chromBPNet** (deep learning models for regulatory prediction)
- **QTL mapping** workflows
- **Figure generation** for manuscript

## Repository Structure

```
cd4t-qtl-map/
├── preprocessing/          # Data preprocessing pipelines
│   ├── ATAC-seq/          # Bulk ATAC-seq processing (Nextflow)
│   └── scRNA-seq/         # Single-cell RNA-seq processing (Nextflow + STARsolo)
├── qtl_mapping/           # QTL analysis workflows
├── chromBPNet_training/   # chromBPNet model training and prediction
├── cd4_qtl_paper_figures/ # Scripts for manuscript figure generation
└── LICENSE
```

## Key Features

### 🧬 Multi-Omic Data Processing
- **ATAC-seq**: Allele-aware alignment (STAR-WASP), peak calling (MACS3), comprehensive QC

    See [`preprocessing/ATAC-seq/README.md`](preprocessing/ATAC-seq/README.md) for:
    - Sample-level preprocessing (FASTQs → peaks + QC)
    - Merged library peak calling and quantification
    - Detailed pipeline descriptions and configurations
- **scRNA-seq**: STARsolo with WASP filtering, demultiplexing, quality control

    See [`preprocessing/scRNA-seq/README.md`](preprocessing/scRNA-seq/README.md) for:
    - Two-part Nextflow workflow (STARsolo + WASP filtering)
    - Matrix generation and demultiplexing
    - Allele-aware processing steps
- **Genotype-aware analysis**: WASP integration for allele-specific mapping

### 📊 Pan CD4+ T cell QTL Mapping
- Integration of chromatin accessibility and gene expression QTLs
- Analysis of immune-related GWAS loci
- Multi-modal regulatory architecture inference

### 🤖 Predictive Modeling
- chromBPNet training for regulatory element prediction in CD4+ T cells
- Model interpretation and variant effect prediction

### 🤖 GWAS colocalization

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- ATAC-seq pipelines adapted from [nf-core/atacseq](https://github.com/nf-core/atacseq)
- scRNA-seq demultiplexing adapted from [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)
- WASP allele-specific alignment from the WASP toolkit
- CellRegMap code and visualization from https://github.com/annacuomo/CellRegMap_analyses 
- MOFA training and visualization adapted from https://biofam.github.io/MOFA2/
- Trackplot scripts and chromBPNet visualizations adapted from https://github.com/GreenleafLab/HDMA, https://github.com/kundajelab/chrombpnet-figures
- Chrombpnet model training pipiline adapted from https://github.com/kundajelab/chrombpnet and children repositories


## Contact

For questions or issues, please open an issue on this repository or contact the authors.

