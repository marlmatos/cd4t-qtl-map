#!/usr/bin/env nextflow

/*
 * ATAC-seq preprocessing pipeline (Part 1)
 * Author: nf-core/atacseq
 * Modified by: Marliette Matos
 *
 * Modification:
 *   - Enables allele-specific alignment using STAR-WASP
 *
 * Date: 2024-08-10
 *
 * Description:
 *   Bulk ATAC-seq preprocessing from FASTQs -> aligned BAMs -> QC -> peak calling per sample.
 *   (Optional) consensus peak generation is included but currently commented out in the workflow.
 */

// -----------------------------------------------------------------------------
// Parameters (paths are currently hard-coded; portability can be addressed later)
// -----------------------------------------------------------------------------
params {
    // Genome resources (see 001_download_hg38.sh / 002_prepare_fasta.sh)
    fasta      = '/gcjm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa'
    fai        = '/gcjm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa.fai'
    gtf        = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf'
    STAR_index = '/gcjm/resources/genome/hg38_gencode_PRI_align'

    bed        = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.bed'
    gene_TSS   = '/gcjm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.tss.bed'
    chromsizes = '/gcjm/resources/genome/hg38.p14.chrom.sizes.fmtd'
    blacklist  = '/gcjm/resources/genome/hg38_gencode_raw/hg38-blacklist.v3.bed'

    // Inputs / outputs
    results    = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/atac_preprocessing'
    data       = '/gcgl/mmatos/cd4_aging_project/data/ATAC-seq/Project_LAP_15324_B01_NAN_Lane.2023-06-06/Sample_*/fastq/*.R{1,2}.fastq.gz'
    datatest   = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/data/cd4_project_atac_fastqs/Project_LAP_15324_B01_NAN_Lane.2023-06-06/Sample_T2612/fastq/*.fastq.gz'

    // Trimming / WASP
    adapter    = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/assets/NexteraPE-PE.fa'
    varVCF     = '/gcjm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4p_WGS_pass_only.snps.maf1.WASP.vcf'  // created with 001_preparing_wasp_vcf.sh

    // Custom scripts/assets
    bam_filter                 = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/bamtools_filtering_script.json'
    bampe_rm_orphan            = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/nfcore_remove_orphans.py'
    macs2_merged_expand_script = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/macs2_merged_expand.py'
    plot_peak_intersect_script = '/gcjm/ATAC-seq_analysis/ATACseq_nf-preprocessing/scripts/bin/plot_peak_intersect.r'

    // Other
    macs_gsize = 'hs'
}

// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------
process FastQC_pretrim {
    tag "$meta.sample_ID"
    publishDir "${params.results}/fastQC/raw_fastqs", mode: 'copy', pattern: "*"

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.{zip,html}', emit: raw_fastqc

    script:
    """
    echo "Performing pre-alignment QC with FastQC"
    fastqc $reads
    """
}

process trim_adapters {
    tag "$meta.sample_ID"
    publishDir "${params.results}/trimgalore", mode: 'copy', pattern: "*"

    input:
    tuple val(meta), path(reads)
    path adapter

    output:
    tuple val(meta), path('*.R{1,2}_val_{1,2}.fq.gz'), emit: trimmed_fastqs

    script:
    """
    echo "Running adapter trimming (trim_galore)"

    trim_galore --fastqc \\
        --basename ${meta.sample_ID} \\
        --paired ${reads} \\
        --adapter file:${adapter} \\
        -O 10 \\
        --output_dir .

    # Rename output files to match expected pattern
    mv "${meta.sample_ID}_val_1.fq.gz" "${meta.sample_ID}.R1_val_1.fq.gz"
    mv "${meta.sample_ID}_val_2.fq.gz" "${meta.sample_ID}.R2_val_2.fq.gz"
    """
}

process FastQC_postrim {
    tag "$meta.sample_ID"
    publishDir "${params.results}/MultiQC/post_trimgalore", mode: 'copy', pattern: "*"

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.html'
    path '*_fastqc.zip'

    script:
    """
    echo "Performing post-trimming QC with FastQC"
    fastqc $reads
    """
}

/*
 * Aligns short sequencing reads with allele-specific resolution (STAR-WASP)
 */

process Genome_Generate {
    publishDir "${params.STAR_index}", mode: 'copy'

    input:
    path fasta
    path fai
    path gtf
    path STAR_index

    output:
    path 'STAR'

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        --genomeDir STAR \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang 100 \\
        --genomeSAindexNbases 12
    """
}

process STAR_align {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/aligned_bams", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(meta), path(reads)
    path index_ch
    path varVCF

    output:
    tuple val(meta), path('*.bam')

    script:
    """
    STAR --runMode alignReads \\
        --genomeDir ${index_ch} \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat \\
        --alignIntronMax 1 \\
        --runThreadN ${task.cpus} \\
        --alignEndsType EndToEnd \\
        --outSAMtype BAM SortedByCoordinate \\
        --runRNGseed 0 \\
        --waspOutputMode SAMtag \\
        --varVCFfile ${varVCF} \\
        --outSAMattributes NH HI AS NM MD vA vG vW \\
        --outFileNamePrefix ${meta.sample_ID}
    """
}

/*
 * Alignment-level QC
 */

process CollectMultipleMetrics {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/Picard_Metrics", mode: 'copy', pattern: "*_metrics"
    publishDir "${params.results}/STAR/Picard_Metrics/plots", mode: 'copy', pattern: "*.pdf"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai

    output:
    tuple val(meta), path("*_metrics"), emit: metrics
    tuple val(meta), path("*.pdf"),     emit: pdf

    script:
    """
    picard \\
        -Xmx1024m \\
        CollectMultipleMetrics \\
        --INPUT $bam \\
        --OUTPUT ${meta.sample_ID}.CollectMultipleMetrics \\
        -R $fasta
    """
}

process Mark_Duplicates {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/Picard_Metrics", mode: 'copy', pattern: "*_metrics.txt"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.MarkDuplicates.bam"), emit: bam
    tuple val(meta), path("*_metrics.txt")

    script:
    """
    picard MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.sample_ID}.MarkDuplicates.bam \\
        --ASSUME_SORTED true \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --REMOVE_DUPLICATES false \\
        --VALIDATION_STRINGENCY LENIENT \\
        --METRICS_FILE ${meta.sample_ID}.marked_dup_metrics.txt
    """
}

process Samtools_Flagstat {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/samtools_stats/pre-filtering", mode: 'copy', pattern: "*.{txt,stats,idxstats}"

    input:
    tuple val(meta), path(bam)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.flagstat.txt")
    tuple val(meta), path("*.stats")
    tuple val(meta), path("*.idxstats")
    tuple val(meta), path("*.bai"), emit: index

    script:
    """
    samtools flagstat ${bam} > ${meta.sample_ID}.bam.flagstat.txt
    samtools stats --threads 1 --reference $fasta ${bam} > ${meta.sample_ID}.sorted.bam.stats
    samtools idxstats --threads 1 ${bam} > ${meta.sample_ID}.sorted.bam.idxstats
    samtools index -@ 1 ${bam}
    """
}

/*
 * Alignment filtering
 */

process Removing_Blacklisted_Regions {
    tag "$meta.sample_ID"

    input:
    tuple val(meta), path(bam), path(bai)
    path blacklist

    output:
    tuple val(meta), path("*.bam")

    script:
    """
    bedtools intersect -v -abam ${bam} -b ${blacklist} > ${meta.sample_ID}.Aligned.sortedbyCoord.MarkDuplicates.blacklst.bam
    """
}

process Post_Alignment_Filtering {
    tag "$meta.sample_ID"

    input:
    tuple val(meta), path(bam)
    path bam_filter

    output:
    tuple val(meta), path("*.bfilt.sorted.bam"), emit: bams

    script:
    """
    samtools index -b ${bam}
    samtools view ${bam} -e 'rname != "chrM"' -o ${meta.sample_ID}.MT.bam

    samtools view -b \\
        -F 0x004 -F 0x0008 \\
        -f 0x001 \\
        -F 0x100 -F 0x0400 -F 0x0200 \\
        -q 1 ${meta.sample_ID}.MT.bam > ${meta.sample_ID}.MT.filt.bam

    samtools view -e '![vW] || [vW]==1' ${meta.sample_ID}.MT.filt.bam -o ${meta.sample_ID}.MT.filt2.bam

    samtools view -h ${meta.sample_ID}.MT.filt2.bam \\
        | awk '{if(\$3 != "chrY" && \$3 != "chrX"){print \$0}}' \\
        | samtools view -Shb - > ${meta.sample_ID}.MT.filt3.bam

    bamtools filter \\
        -in ${meta.sample_ID}.MT.filt3.bam \\
        -out ${meta.sample_ID}.bfilt.bam \\
        -script $bam_filter

    samtools sort -@ 9 -o ${meta.sample_ID}.bfilt.sorted.bam ${meta.sample_ID}.bfilt.bam
    """
}

process Remove_Orphans {
    tag "$meta.sample_ID"

    input:
    tuple val(meta), path(bam)
    path bampe_rm_orphan

    output:
    tuple val(meta), path('*.bfilt.sorted.orpham.sorted.bam'), emit: filt_bams

    script:
    """
    nfcore_remove_orphans.py $bam ${meta.sample_ID}.bfilt.sorted.orpham.bam --only_fr_pairs
    """
}

process Bam_Index {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/filtered_bams", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.results}/STAR/samtools_stats/post-filtering", mode: 'copy', pattern: "*.{txt,stats,idxstats}"

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: sorted_bams
    tuple val(meta), path("*.flagstat.txt"), emit: flagstat
    tuple val(meta), path("*.stats")
    tuple val(meta), path("*.idxstats")

    script:
    """
    samtools flagstat ${bam} > ${meta.sample_ID}.bfilt.sorted.orpham.bam.flagstat.txt
    samtools stats --threads 1 --reference $fasta ${bam} > ${meta.sample_ID}.bfilt.sorted.orpham.bam.stats
    samtools idxstats --threads 1 ${bam} > ${meta.sample_ID}.bfilt.sorted.orpham.bam.idxstats

    samtools sort -@ 9 -o ${meta.sample_ID}.bfilt.sorted.orpham.sorted.bam ${bam}
    samtools index ${meta.sample_ID}.bfilt.sorted.orpham.sorted.bam
    """
}

/*
 * ATAC shift + indexing
 */

process Shift_Alignment {
    tag "$meta.sample_ID"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.shifted_filtered.bam"), emit: shifted_bam

    script:
    """
    alignmentSieve -b $bam \\
        --ATACshift \\
        -o ${meta.sample_ID}.shifted_filtered.bam
    """
}

process Index_Shifted_Bam {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/shifted_bams", mode: 'copy'

    input:
    tuple val(meta), path(shifted_bam)

    output:
    tuple val(meta), path("*.shifted.sorted.bam"), path("*.shifted.sorted.bam.bai"), emit: bams

    script:
    """
    samtools sort -@ 9 -o ${meta.sample_ID}.shifted.sorted.bam ${shifted_bam}
    samtools index -b ${meta.sample_ID}.shifted.sorted.bam
    """
}

/*
 * Fingerprint QC
 */

process Fingerprint_profile {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/deeptools", mode: 'copy', pattern: "*.txt"
    publishDir "${params.results}/STAR/deeptools/plots", mode: 'copy', pattern: "*.pdf"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.pdf")
    tuple val(meta), path("*.raw.txt")
    tuple val(meta), path("*.qcmetrics.txt")

    script:
    """
    plotFingerprint \\
        --skipZeros --numberOfSamples 500000 --labels ${meta.sample_ID} \\
        --bamfiles ${bam} \\
        --plotFile ${meta.sample_ID}.plotFingerprint.pdf \\
        --outRawCounts ${meta.sample_ID}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${meta.sample_ID}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors 12
    """
}

/*
 * Peak calling + signal tracks
 */

process MACS3_Peakcall {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/MACS3/peaks",   mode: 'copy', pattern: "*.{narrowPeak,xls}"
    publishDir "${params.results}/STAR/MACS3/pileups", mode: 'copy', pattern: "*.bdg"

    input:
    tuple val(meta), path(shifted_bam)
    val macs_gsize

    output:
    tuple val(meta), path("*_treat_pileup.bdg"),   emit: bedgraph
    tuple val(meta), path("*_control_lambda.bdg"), emit: lambda_bedgraph
    tuple val(meta), path("*.xls"),                emit: peaks_xls
    tuple val(meta), path("*_peaks.narrowPeak"),   emit: narrow_peaks

    script:
    """
    macs3 callpeak \\
        -t $shifted_bam \\
        -f BAMPE \\
        -n $meta.sample_ID \\
        -g $macs_gsize \\
        --nomodel \\
        -B \\
        -q 0.01 \\
        --nolambda \\
        --seed 762873 \\
        --keep-dup all \\
        --SPMR
    """
}

process BedGraphToBigWig {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/bigwig", mode: 'copy', pattern: "*.bigWig"

    input:
    tuple val(meta), path(bedgraph)
    path chromsizes

    output:
    tuple val(meta), path("*.bigWig"), emit: bigWig

    script:
    """
    bedGraphToBigWig \\
        $bedgraph \\
        $chromsizes \\
        ${meta.sample_ID}.bigWig
    """
}

process Deeptools_QC_Scaled_Regions {
    tag "$meta.sample_ID"
    publishDir "${params.results}/STAR/deeptools",       mode: 'copy', pattern: "*.{gz,tab}"
    publishDir "${params.results}/STAR/deeptools/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple val(meta), path(bigWig)
    path bed

    output:
    tuple val(meta), path("*.mat.gz")
    tuple val(meta), path("*.tab")
    tuple val(meta), path("*.png")

    script:
    """
    computeMatrix scale-regions \\
        -R $bed \\
        -S $bigWig \\
        -o ${meta.sample_ID}.scaleregions.computeMatrix.mat.gz \\
        --outFileNameMatrix ${meta.sample_ID}.scale_regions.computeMatrix.vals.mat.tab \\
        -b 3000 -a 3000 \\
        --regionBodyLength 1000 \\
        --missingDataAsZero \\
        --skipZeros \\
        --smartLabels \\
        --numberOfProcessors $task.cpus

    plotProfile \\
        --matrixFile ${meta.sample_ID}.scaleregions.computeMatrix.mat.gz \\
        --outFileName ${meta.sample_ID}.scaleregions.plotProfile.png \\
        --outFileNameData ${meta.sample_ID}.scaleregions.plotProfile.tab
    """
}

process Deeptools_QC_Reference {
    tag "merged_library"
    publishDir "${params.results}/STAR/deeptools",       mode: 'copy', pattern: "*.gz"
    publishDir "${params.results}/STAR/deeptools/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple val(meta), path(bigWig)
    path gene_TSS

    output:
    tuple val(meta), path("*.mat.gz")
    tuple val(meta), path("*.png")

    script:
    """
    computeMatrix reference-point \\
        --missingDataAsZero --skipZeros -b 3000 -a 3000 --smartLabels \\
        -R $gene_TSS \\
        -S $bigWig \\
        -o ${meta.sample_ID}.referencepoint.computeMatrix.mat.gz \\
        --numberOfProcessors $task.cpus

    plotHeatmap \\
        --matrixFile ${meta.sample_ID}.referencepoint.computeMatrix.mat.gz \\
        --outFileName ${meta.sample_ID}.referencepoint.plotHeatmap.png
    """
}

/*
 * Consensus peaks (optional; currently commented out in workflow)
 */

process Consensus_peaks {
    publishDir "${params.results}/STAR/MACS3/Concensus", mode: 'copy'

    input:
    tuple val(meta), path(narrow_peaks)
    path macs2_merged_expand_script
    path plot_peak_intersect_script

    output:
    tuple val(meta), path("*.bed"),           emit: bed
    tuple val(meta), path("*.saf"),           emit: saf
    tuple val(meta), path("*.pdf"),           emit: pdf
    tuple val(meta), path("*.boolean.txt"),   emit: boolean_txt
    tuple val(meta), path("*.intersect.txt"), emit: intersect_txt

    script:
    def args         = "--min_replicates 1"
    def prefix       = "cd4_atac_consensus_peaks"
    def peak_type    = "narrowPeak"
    def mergecols    = (2..10).join(',')
    // 9 columns collapsed across cols 2..10
    def collapsecols = (['collapse'] * 9).join(',')
    def expandparam  = "--is_narrow_peak"
    def fdr          = "--qvalue-threshold 0.01"

    """
    sort -T '.' -k1,1 -k2,2n ${narrow_peaks.collect{ it.toString() }.sort().join(' ')} \\
        | mergeBed -c $mergecols -o $collapsecols > ${prefix}.txt

    python ${macs2_merged_expand_script} \\
        ${prefix}.txt \\
        ${narrow_peaks.collect{ it.toString() }.sort().join(',').replaceAll("_peaks.${peak_type}","")} \\
        ${prefix}.boolean.txt \\
        $args \\
        $expandparam \\
        $fdr

    awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 { print \$1, \$2, \$3, \$4, "0", "+" }' ${prefix}.boolean.txt > ${prefix}.bed

    echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > ${prefix}.saf
    awk -v FS='\\t' -v OFS='\\t' 'FNR > 1 { print \$4, \$1, \$2, \$3, "+" }' ${prefix}.boolean.txt >> ${prefix}.saf

    Rscript ${plot_peak_intersect_script} ${prefix}.boolean.txt ${prefix}.pdf
    """
}

// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------
workflow {

    trimmed_reads_ch = Channel
        .fromFilePairs(params.data)
        .map { id, reads ->
            def meta = [
                ['sample_ID', 'i7-i5_indexes', 'flowcell_ID', 'lane', 'replicate'],
                id.tokenize("_")
            ].transpose().collectEntries()
            [meta, reads]
        }
        .set { fastqs }

    FastQC_pretrim(fastqs)

    trim_adapters(fastqs, params.adapter)
    FastQC_postrim(trim_adapters.out.trimmed_fastqs)

    Genome_Generate(params.fasta, params.fai, params.gtf, params.STAR_index)

    STAR_align(trim_adapters.out.trimmed_fastqs, Genome_Generate.out, params.varVCF)

    CollectMultipleMetrics(STAR_align.out, params.fasta, params.fai)
    Mark_Duplicates(STAR_align.out, params.fasta)

    Samtools_Flagstat(Mark_Duplicates.out[0], params.fasta, params.fai)

    // Post-alignment filtering
    Mark_Duplicates.out.bam
        .combine(Samtools_Flagstat.out.index, by: 0)
        .view()
        .set { blacklist_input_ch }

    Removing_Blacklisted_Regions(blacklist_input_ch, params.blacklist)

    Post_Alignment_Filtering(Removing_Blacklisted_Regions.out, params.bam_filter)

    Remove_Orphans(Post_Alignment_Filtering.out.bams, params.bampe_rm_orphan)

    Bam_Index(Remove_Orphans.out.filt_bams, params.fasta)

    // Shift alignments before peak calling so peaks match signal tracks
    Shift_Alignment(Bam_Index.out.sorted_bams)

    // QC
    Index_Shifted_Bam(Shift_Alignment.out.shifted_bam)
    Fingerprint_profile(Index_Shifted_Bam.out.bams)

    Shift_Alignment.out.shifted_bam
        .combine(Bam_Index.out.flagstat, by: 0)
        .view()
        .set { gcov_input_ch }

    // Peak calling
    MACS3_Peakcall(Shift_Alignment.out, params.macs_gsize)

    BedGraphToBigWig(MACS3_Peakcall.out.bedgraph, params.chromsizes)

    Deeptools_QC_Scaled_Regions(BedGraphToBigWig.out.bigWig, params.bed)
    Deeptools_QC_Reference(BedGraphToBigWig.out.bigWig, params.gene_TSS)

    MACS3_Peakcall.out.narrow_peaks
        .combine(gcov_input_ch, by: 0)
        .view()
        .set { frip_input_ch }

    // Optional downstream steps (currently disabled)
    // Calculate_FRIP_1(frip_input_ch)
    //
    // MACS3_Peakcall.out.narrow_peaks
    //     .collect { it[1] }
    //     .filter { it.size() > 1 }
    //     .map { peaks -> [ [ id: 'consensus_peaks' ], peaks ] }
    //     .view()
    //     .set { ch_consensus_peaks }
    //
    // Consensus_peaks(ch_consensus_peaks, params.macs2_merged_expand_script, params.plot_peak_intersect_script)
    // Annotate_peaks_1(Consensus_peaks.out.bed, params.fasta, params.gtf)
    // Shift_Alignment.out.shifted_bam.combine(Consensus_peaks.out.saf, by: 0).view().set { Feature_Counts_input_ch }
    // Feature_Counts(Feature_Counts_input_ch)
}

