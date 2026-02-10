# CD4 ATAC-seq preprocessing (Nextflow)

Nextflow pipelines and helper scripts to preprocess bulk ATAC-seq data from FASTQs through alignment (STAR-WASP), QC, filtering, shifting, peak calling (MACS3), and downstream QC tracks/metrics. This repository contains the code used for figure generation and analyses associated with:

**Integrating multi-omic QTLs and predictive models reveals regulatory architectures at immune related GWAS loci in CD4+ T cells.**

**Please note:** most scripts used in this repo are modifications from the [nextflow/atac-seq](https://github.com/nf-core/atacseq)
---

## What’s in here

### Pipeline 1: Sample-level preprocessing (FASTQs → peaks + QC)
Main steps:
- **FastQC** (pre/post trim)
- **Trim Galore** adapter trimming
- **STAR** alignment with **WASP** tags for allele-specific alignment
- **Picard** metrics + MarkDuplicates
- **Samtools** stats/flagstat/index
- **Blacklist filtering** (bedtools intersect -v)
- **Post-alignment filtering** (bamtools JSON script; removes chrM, MAPQ filters, vW/WASP, etc.)
- **Orphan removal** (nf-core helper)
- **Tn5 shift** (deepTools `alignmentSieve --ATACshift`)
- **Fingerprint QC** (deepTools `plotFingerprint`)
- **Peak calling** (MACS3)
- **BigWig generation** (bedGraphToBigWig)
- **DeepTools QC** around gene bodies and TSS

### Pipeline 2: Merged-library peak calling and quantification (shifted BAMs → merged peaks)
Main steps:
- Merge shifted BAMs into a **merged library**
- Sort/index merged BAM
- MACS3 peak calling (optionally call summits)
- Pad/standardize peaks around summits (custom script)
- FRiP against processed peaks
- BigWig + deepTools QC for merged library
- FeatureCounts quantification over processed peaks

---

<h2>Repository structure</h2>
<pre><code>.
├── assets/
├── bin/
├── conf/
├── main.groovy
├── main_merged_library.groovy
├── nextflow.config
├── run_nextflow_atac_pt1.sh
└── run_nextflow_atac_pt2.sh
</code></pre>

## Description
- **`main.groovy`**  
  Pipeline Part 1: sample-level preprocessing (FASTQs → aligned/filtered/shifted BAMs + MACS3 peaks + QC).

- **`main_merged_library.groovy`**  
  Pipeline Part 2: merged-library workflow (shifted BAMs → merged BAM → merged peak calling + processed peaks + FRiP + featureCounts + QC).

- **`nextflow.config`**  
  Cluster execution config (e.g., Slurm), resources, retry strategy, Singularity container definitions, and any profile settings.

- **`conf/`**  
  Optional: extra config files (cluster profiles, resource overrides, etc.). If you have multiple environments, put them here.

- **`assets/`**  
  Static inputs referenced by the pipelines (e.g., adapter fasta, headers, templates).

- **`bin/`**  
  Small helper scripts invoked by pipeline processes (Python/R scripts, JSON filter specs).

- **`run_nextflow_atac_pt1.sh` / `run_nextflow_atac_pt2.sh`**  
  Convenience wrappers to run each part with the right profile/config and common parameters.
