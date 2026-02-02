#!/usr/bin/env nextflow

/*
 * Bulk ATAC-seq preprocessing (Merged-library peak calling)
 * FASTQs -> aligned/shifted BAMs (produced elsewhere) -> merged BAM -> MACS3 peaks -> QC -> FRiP -> featureCounts -> Homer annotation
 *
 * NOTE: This script currently consumes:
 *   - shifted BAMs (params.shifted_bams)
 *   - corresponding flagstat files (params.shifted_bams_flagstat)
 * and performs merged-library peak calling + downstream QC/quantification.
 */

// -----------------------------------------------------------------------------
// Parameters (resources, inputs, outputs)
// -----------------------------------------------------------------------------

// Resources (see: 001_download_hg38.sh and 002_prepare_fasta.sh)
params.fasta      = '/gcjm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa'
params.fai        = '/gcjm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa.fai'
params.gtf        = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf'
params.STAR_index = '/gcjm/resources/genome/hg38_gencode_PRI_align'

params.bed        = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.bed'
params.gene_TSS   = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.tss.bed'
params.chromsizes = '/gcjm/resources/genome/hg38.p14.chrom.sizes.fmtd'
params.blacklist  = '/gcjm/resources/genome/hg38_gencode_raw/hg38-blacklist.v3.bed'

// Inputs / outputs
params.results   = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing'
params.data      = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/data/cd4_project_atac_fastqs/Project_LAP_15324_B01_NAN_Lane.2023-06-06/Sample_*/fastq/*.R{1,2}.fastq.gz'
params.datatest  = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/data/cd4_project_atac_fastqs/Project_LAP_15324_B01_NAN_Lane.2023-06-06/Sample_T2612/fastq/*.fastq.gz'

// Assets / helper scripts
params.adapter                   = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/assets/NexteraPE-PE.fa'
params.varVCF                    = '/gcjm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4p_WGS_pass_only.snps.maf1.WASP.vcf' // created with 001_preparing_wasp_vcf.sh
params.bam_filter                = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/bamtools_filtering_script.json'
params.bampe_rm_orphan           = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/nfcore_remove_orphans.py'
params.macs2_merged_expand_script = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/macs2_merged_expand.py'
params.plot_peak_intersect_script = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/plot_peak_intersect.R'
params.peaks_script              = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/unique_peaks.py'
params.custompeaks_script        = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/customPeakset.r'
params.igv_script                = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/igv_files_to_session.py'
params.pad_peaks_script          = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/pad_peaks.py'

params.rlibrary          = '/gcjm/R/x86_64-pc-linux-gnu-library'
params.hommer_header     = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/assets/merged_library_peak_annotation_header.txt'
params.hommer_plot_script = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/plot_homer_annotatepeaks.r'

// Precomputed outputs from upstream pipeline
params.shifted_bams          = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/shifted_bams/*.bam'
params.shifted_bams_flagstat = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/samtools_stats/post-filtering/*.bfilt.orpham.bam.flagstat.txt'
params.bigwig                = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing/STAR/bigwig/*.bigWig'

// Peak processing params
params.extension       = 300
params.score_threshold = 0.05

// Other
params.macs_gsize = 'hs'


// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------

/*
 * Merge shifted BAMs into a single mega-library BAM for joint peak calling.
 */
process MergeBams {
    tag "merged_library"

    input:
    path(bam)

    output:
    path("*.list"),    emit: merge_list
    path("*.out.bam"), emit: merged_bam

    script:
    """
    echo "${bam.join('\n')}" > bam.list
    samtools merge --threads ${task.cpus} -f -o merged.out.bam -b bam.list
    """
}

/*
 * Sort + index merged BAM; compute merged flagstat.
 */
process Sort_Index {
    tag "merged_library"
    publishDir "${params.results}/merged_library/bam", mode: 'copy', pattern: "*.{bam,bai,txt}"

    input:
    path(merged_bam_ch)

    output:
    path("*.sorted.bam"),     emit: sorted_bam
    path("*.sorted.bam.bai"), emit: sorted_bai
    path("*.txt"),            emit: merged_flagstat

    script:
    def prefix = "cd4_atac"
    """
    samtools flagstat ${merged_bam_ch} > ${prefix}.merged.out.flagstat.txt
    samtools sort -@ ${task.cpus} -o ${prefix}.merged.out.sorted.bam merged.out.bam
    samtools index ${prefix}.merged.out.sorted.bam
    """
}

/*
 * Convert merged BAM to BED (sometimes useful for debugging / downstream scripts).
 */
process BAM_to_BED {
    tag "merged_library"

    input:
    path(sorted_bam)

    output:
    path("*.bed"), emit: bed

    script:
    def prefix = "cd4_atac"
    """
    bedtools bamtobed -i ${sorted_bam} > ${prefix}_merged_lib.bed
    """
}

/*
 * Peak calling on merged library (MACS3; BAMPE).
 */
process MACS2_Peakcall_1 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/peaks_102024",   mode: 'copy', pattern: "*.{narrowPeak,xls,txt}"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/pileups_102024", mode: 'copy', pattern: "*.bdg"

    input:
    path(sorted_bam)
    val(macs_gsize)

    output:
    path("*_treat_pileup.bdg"),    emit: bedgraph
    path("*_control_lambda.bdg"),  emit: lambda_bedgraph
    path("*.xls"),                 emit: peaks_xls
    path("*_peaks.narrowPeak"),    emit: narrow_peaks
    path("*_summits.bed"),         emit: cutoff_analysis

    script:
    def prefix = "cd4_atac"
    """
    macs3 callpeak \\
      -t ${sorted_bam} \\
      -f BAMPE \\
      -n ${prefix} \\
      -g ${macs_gsize} \\
      --nomodel \\
      -B \\
      --min-length 180 \\
      --max-gap 1000 \\
      -q 0.01 \\
      --nolambda \\
      --seed 762873 \\
      --keep-dup all \\
      --SPMR
    """
}

/*
 * Peak calling with summits.
 */
process MACS2_Peakcall_2 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/peaks_102024",   mode: 'copy', pattern: "*.{narrowPeak,xls,txt}"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/pileups_102024", mode: 'copy', pattern: "*.bdg"

    input:
    path(sorted_bam)
    val(macs_gsize)

    output:
    path("*summits_treat_pileup.bdg"), emit: bedgraph
    path("*_control_lambda.bdg"),      emit: lambda_bedgraph
    path("*.xls"),                     emit: peaks_xls
    path("*_peaks.narrowPeak"),        emit: narrow_peaks
    path("*_summits.bed"),             emit: cutoff_analysis

    script:
    def prefix = "cd4_atac_summits"
    """
    macs3 callpeak \\
      -t ${sorted_bam} \\
      -f BAMPE \\
      -n ${prefix} \\
      -g ${macs_gsize} \\
      --nomodel \\
      -B \\
      --min-length 180 \\
      --max-gap 1000 \\
      -q 0.01 \\
      --nolambda \\
      --seed 762873 \\
      --keep-dup all \\
      --SPMR \\
      --call-summits
    """
}

/*
 * Convert bedGraph (MACS pileup) -> bigWig for visualization and deepTools QC.
 */
process BedGraphToBigWig_2 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/bigwig", mode: 'copy', pattern: "*.bigWig"

    input:
    path(bedgraph)
    path(chromsizes)

    output:
    path("*.bigWig"), emit: bigWig

    script:
    def prefix = "cd4_atac"
    """
    bedGraphToBigWig \\
      ${bedgraph} \\
      ${chromsizes} \\
      ${prefix}.bigWig
    """
}

/*
 * deepTools: computeMatrix + plotProfile over annotation BED.
 */
process Deeptools_QC_Scaled_Regions_2 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/deeptools",       mode: 'copy', pattern: "*.{gz,tab}"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/deeptools/plots", mode: 'copy', pattern: "*.png"

    input:
    path(bigWig)
    path(bed)

    output:
    path("*.mat.gz")
    path("*.tab")
    path("*.png")

    script:
    def prefix = "cd4_atac"
    """
    computeMatrix scale-regions \\
      -R ${bed} \\
      -S ${bigWig} \\
      -o ${prefix}.scaleregions.computeMatrix.mat.gz \\
      --outFileNameMatrix ${prefix}.scale_regions.computeMatrix.vals.mat.tab \\
      -b 3000 \\
      -a 3000 \\
      --regionBodyLength 1000 \\
      --missingDataAsZero \\
      --skipZeros \\
      --smartLabels \\
      --numberOfProcessors ${task.cpus}

    plotProfile \\
      --matrixFile ${prefix}.scaleregions.computeMatrix.mat.gz \\
      --outFileName ${prefix}.scaleregions.plotProfile.png \\
      --outFileNameData ${prefix}.scaleregions.plotProfile.tab
    """
}

/*
 * deepTools: reference-point heatmap around TSS.
 */
process Deeptools_QC_Reference_2 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/deeptools",       mode: 'copy', pattern: "*.gz"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/deeptools/plots", mode: 'copy', pattern: "*.png"

    input:
    path(bigWig)
    path(gene_TSS)

    output:
    path("*.mat.gz")
    path("*.png")

    script:
    def prefix = "cd4_atac"
    """
    computeMatrix reference-point --missingDataAsZero --skipZeros -b 3000 -a 3000 --smartLabels \\
      -R ${gene_TSS} \\
      -S ${bigWig} \\
      -o ${prefix}.referencepoint.computeMatrix.mat.gz \\
      --numberOfProcessors ${task.cpus}

    plotHeatmap \\
      --matrixFile ${prefix}.referencepoint.computeMatrix.mat.gz \\
      --outFileName ${prefix}.referencepoint.plotHeatmap.png
    """
}

/*
 * Post-process MACS summits: pad peaks, filter by score, respect chrom sizes.
 */
process PROCESS_MACS2_SUMMITS {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/peaks_102024", mode: 'copy', pattern: "*.bed"
    conda "/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/python_env.yml"

    input:
    path(narrow_peaks)
    path(chromsizes)
    val(extension)
    val(score_threshold)
    path(pad_peaks_script)

    output:
    path("padded_summits_peaks.bed"), emit: processed_peaks

    script:
    """
    python pad_peaks.py \\
      --input ${narrow_peaks} \\
      --output padded_summits_peaks.bed \\
      --chrom_sizes ${chromsizes} \\
      --extension ${extension} \\
      --score_threshold ${score_threshold} \\
      --verbose
    """
}

/*
 * FRiP per sample using processed peaks.
 */
process Calculate_FRIP_2 {
    tag "${meta.sample}"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/FRIP/", mode: 'copy', pattern: "*.FRiP.txt"

    input:
    tuple val(meta), path(bam), path(flagstat)
    path(processed_peaks)

    output:
    path("*.FRiP.txt"), emit: FRiP_narrow_peaks

    script:
    """
    READS_IN_PEAKS=\$(intersectBed -a ${bam} -b ${processed_peaks} -bed -c -f 0.20 -sorted | awk -F '\\t' '{sum += \$NF} END {print sum}')
    grep 'mapped (' ${flagstat} | grep -v "primary" | awk -v a="\$READS_IN_PEAKS" -v OFS='\\t' '{print "${meta.sample}", a/\$1}' > ${meta.sample}.narrow_peaks.FRiP.txt
    """
}

/*
 * Homer annotation of processed peaks.
 */
process Annotate_peaks_2 {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/annotate_peaks", mode: 'copy', pattern: "*.{narrowPeak,xls,txt}"

    input:
    path(processed_peaks)
    path(fasta)
    path(gtf)

    output:
    path("*annotatePeaks.txt"), emit: txt
    path("*annStats.txt"),      emit: stats, optional: true

    script:
    def prefix = "cd4_atac"
    """
    annotatePeaks.pl \\
      ${processed_peaks} \\
      ${fasta} \\
      -gid \\
      -gtf ${gtf} \\
      -cpu ${task.cpus} \\
      > ${prefix}_narrowPeaks.annotatePeaks.txt
    """
}

/*
 * featureCounts over processed peaks SAF; counts across all shifted BAMs.
 */
process Feature_Counts_merged {
    tag "merged_library"
    publishDir "${params.results}/merged_library/peak_calling/MACS3/BAMPE/feature_counts", mode: 'copy'

    input:
    path(shifted_bam_ch)
    path(processed_peaks)

    output:
    path("*featureCounts.txt")        , emit: counts
    path("*featureCounts.txt.summary"), emit: summary, optional: true
    path("*narrowPeaks.saf")          , emit: saf

    script:
    def prefix = "cd4_atac"
    """
    echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > ${prefix}_narrowPeaks.saf
    awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 { print \$4, \$1, \$2, \$3, "+" }' ${processed_peaks} >> ${prefix}_narrowPeaks.saf

    featureCounts \\
      -F SAF \\
      -O --fracOverlap 0.2 \\
      -p \\
      -T ${task.cpus} \\
      -a ${prefix}_narrowPeaks.saf \\
      -s 0 \\
      -o ${prefix}_narrowPeaks.featureCounts.txt \\
      ${shifted_bam_ch.join(' ')}
    """
}


// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------

workflow {

    // Collect shifted BAMs (from upstream pipeline outputs)
    Channel.fromPath(params.shifted_bams)
        | collect
        | set { shifted_bam_ch }

    MergeBams(shifted_bam_ch)
    Sort_Index(MergeBams.out.merged_bam)

    // Peak calling on merged library
    BAM_to_BED(Sort_Index.out.sorted_bam)

    MACS2_Peakcall_1(Sort_Index.out.sorted_bam, params.macs_gsize)
    MACS2_Peakcall_2(Sort_Index.out.sorted_bam, params.macs_gsize)

    // Process summits into padded peak set
    PROCESS_MACS2_SUMMITS(
        MACS2_Peakcall_2.out.narrow_peaks,
        params.chromsizes,
        params.extension,
        params.score_threshold,
        params.pad_peaks_script
    )

    // Build a channel of flagstat files with meta = [sample: <prefix>]
    Channel.fromPath(params.shifted_bams_flagstat)
        .map { file ->
            def sample = file.simpleName
            def meta   = [ sample: sample ]
            [ meta, file ]
        }
        .transpose()
        .set { flagstat_ch }

    // Build a channel of shifted BAMs with meta = [sample: <prefix>]
    shifted_bam_ch
        .flatten()
        .map { file ->
            def sample = file.simpleName
            def meta   = [ sample: sample ]
            [ meta, file ]
        }
        .set { processed_shifted_bam_ch }

    // Join BAM + flagstat by sample ID (meta map)
    processed_shifted_bam_ch
        .combine(flagstat_ch, by: 0)
        .view()
        .set { frip_input_ch }

    // FRiP per sample
    Calculate_FRIP_2(frip_input_ch, PROCESS_MACS2_SUMMITS.out.processed_peaks)

    // bigWig + deepTools QC from merged library pileup bedgraph
    BedGraphToBigWig_2(MACS2_Peakcall_2.out.bedgraph, params.chromsizes)
    Deeptools_QC_Scaled_Regions_2(BedGraphToBigWig_2.out.bigWig, params.bed)
    Deeptools_QC_Reference_2(BedGraphToBigWig_2.out.bigWig, params.gene_TSS)

    // Quantify peaks across all shifted BAMs
    Feature_Counts_merged(shifted_bam_ch, PROCESS_MACS2_SUMMITS.out.processed_peaks)

    // Annotate peaks
    Annotate_peaks_2(PROCESS_MACS2_SUMMITS.out.processed_peaks, params.fasta, params.gtf)

    // Optional:
    // bigwig_ch = Channel.fromPath(params.bigwig).collect()
    // Annotate_peaks_QC(Annotate_peaks_2.out.txt, params.hommer_header, params.hommer_plot_script)
    // IGV_session(params.fasta, params.fai, params.igv_script, bigwig_ch, MACS2_Peakcall_2.out.narrow_peaks, MACS2_Peakcall_2.out.bedgraph)
}
