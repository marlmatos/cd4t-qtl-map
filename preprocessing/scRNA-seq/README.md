# scRNA-seq preprocessing workflow (pt1 + pt2)

This directory contains a two-part Nextflow workflow for preprocessing 10x Genomics scRNA-seq libraries using **STARsolo** with **STAR-WASP** allele-aware alignment.

**Important design note (why there are two parts):**
STARsolo can add WASP tags to alignments (e.g., `vW`), but the **default STARsolo count matrix is generated without filtering on `vW`**. In our workflow, downstream analyses require removing reads that fail WASP filtering. Therefore, the pipeline is split into:

- **pt1 (`scRNAseq_preprocessing.groovy`)**: run STARsolo alignment/quantification with WASP tags enabled and generate the per-library outputs (BAM + STARsolo directory).
- **pt2 (`scRNAseq_preprocessing2.groovy`)**: **filter the STARsolo BAMs using the WASP tag (`vW==1` or missing `vW`)**, restrict to autosomes, then **recreate the gene-by-cell count matrix** from the filtered alignments and proceed to downstream steps (e.g., demultiplexing).

#### Please note:
Parts of this pipeline (cluster-to-donor assignment via genotype correlation and related helper scripts)
are adapted from Demuxafy (https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/) . We modified the original code for our study-specific inputs/outputs and
integration into this Nextflow workflow.

In other words, **pt2 is not optional** in the WASP-aware workflow: it is the step that produces the **analysis-ready matrices**.

---

<h2>Directory structure</h2>
<pre><code>scRNAseq/
├── assets/
│   └── Static resources used by the pipeline (e.g., barcode whitelist)
├── bin/
│   └── Helper scripts called by Nextflow processes (Python/R; executable)
├── scRNAseq_preprocessing.groovy
│   └── pt1: STARsolo (+ WASP tags) alignment/quantification outputs
├── scRNAseq_preprocessing2.groovy
│   └── pt2: WASP-filter BAMs → rebuild matrices → downstream steps
├── nextflow.config
│   └── Global Nextflow configuration + profiles + container/module settings
├── nextflow_scrnaseq_STAR.sh
│   └── SLURM wrapper (currently configured to run pt2; edit the `nextflow run`
│       target to switch between pt1/pt2)
├── 001_download_hg38.sh / 002_prepare_fasta.sh
│   └── Reference genome download/prep
├── 001_renaming_fastq.sh / 001_renaming_fastqs.py
│   └── Standardize raw FASTQ filenames before running pt1
├── 001_preparing_wasp_vcf.sh
│   └── Build the WASP VCF used for STAR-WASP allele-aware alignment (main script)
└── 002_preparing_demultiplexing_vcf.sh
    └── Build the multi-sample donor VCF used for demultiplexing/correlation
</code></pre>


---

## What each part produces

### pt1 (`scRNAseq_preprocessing.groovy`) — STARsolo alignment with WASP tags
Typical input: **FASTQs**.

Key outputs (per library under `${outdir}/STAR/<library_id>/`):
- `Aligned.sortedByCoord.out.bam` (contains WASP tags, including `vW`)
- `<library_id>.Solo.out/` (STARsolo output directory; includes *unfiltered* raw matrices)
- `Log.final.out` and other STAR logs

**Note:** The STARsolo matrix produced here is **not** WASP-filtered and is **not used for final downstream analysis** in this workflow.

---

### pt2 (`scRNAseq_preprocessing2.groovy`) — WASP-aware filtering + matrix rebuild (required)
Typical input: **STAR outputs from pt1** (BAM + `Solo.out` directory).

Key steps:
1. **Filter alignments by WASP tag**
   - Keep reads that **passed WASP** (`vW==1`) and reads with **no `vW` tag** (no overlapping variant)
2. Filter to autosomes (remove `chrX`, `chrY`)
3. Convert filtered BAM → SAM
4. **Rebuild `matrix.mtx`, `barcodes.tsv`, `features.tsv`** from filtered alignments
   - Uses STARsolo `barcodes.tsv` / `features.tsv` plus a BAM→matrix reconstruction script
5. Proceed to downstream analysis steps (e.g., souporcell demultiplexing and correlation)

Key outputs (per library):
- filtered BAM + index
- **analysis-ready** rebuilt count matrices from filtered reads
- demultiplexing outputs (if enabled)

---

## How pt1 and pt2 connect

- pt1 writes per-library STAR outputs to:
  - `${outdir}/STAR/<library_id>/`
- pt2 reads those directories via something like:
  - `params.STAR_align_output = "${outdir}/STAR/scGEX-A*"`
  and expects (per library):
  - `${paths}/${sample_ID}.Aligned.sortedByCoord.out.bam`
  - `${paths}/${sample_ID}.Solo.out/Gene/raw/{barcodes.tsv, features.tsv}`

---

## Required files

### Core pipeline files
- `scRNAseq_preprocessing.groovy` (pt1)
- `scRNAseq_preprocessing2.groovy` (pt2)
- `nextflow.config`
- `nextflow_scrnaseq_STAR.sh`

### Helper scripts (in `bin/`)
Used to rebuild matrices / demultiplex / assign clusters:
- `soloCOUNTmatrixFromBam.py` — rebuild the count matrix from WASP-filtered alignments
- `header_matrix.py` — write matrix headers
- `batch_txt_script.r` — batch helper for demultiplexing
- `demulti_corr.r` — genotype–cluster correlation (donor assignment)

### Required reference resources (passed via `params.*`)
- Genome: `fasta`, `fai`, `gtf`, `STAR_index/`
- 10x barcode whitelist
- WASP VCF for STAR alignment
- Donor multi-sample VCF + metadata CSV (for demultiplexing/correlation steps)

