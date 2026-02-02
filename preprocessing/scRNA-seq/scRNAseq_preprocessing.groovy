```nextflow
#!/usr/bin/env nextflow

/*
 * scRNA-seq preprocessing pipeline (Part 1): alignment with STARsolo + WASP and initial count generation
 *
 * Based on: nf-core/scrnaseq (adapted)
 * Modified by: Marliette Matos
 * Date: <2024-04-01>
 *
 * Purpose
 * -------
 * This workflow preprocesses 10x Genomics 3' Gene Expression (v3) scRNA-seq FASTQs through:
 *   (1) STAR genome index generation (if needed)
 *   (2) STARsolo alignment and gene counting with allele-aware alignment enabled via STAR WASP mode
 *   (3) Post-alignment BAM filtering to retain autosomal reads and WASP-passing reads
 *   (4) Rebuilding the raw Gene count matrix (matrix.mtx) from the filtered alignments
 *   (5) EmptyDrops-based cell filtering to produce filtered matrices
 *
 * Major steps
 * -----------
 * 1) Genome_Generate: build STAR index from FASTA + GTF
 * 2) STAR_align: align reads with STARsolo (10x v3), output BAM + Solo.out gene counts; enable WASP tags (vW)
 * 3) (Optional) STAR_multiqc: summarize STAR Log.final.out across libraries
 * 4) Bam_Filter: filter BAM for autosomes and WASP-passing reads (vW==1 or missing vW tag)
 * 5) BAM_to_SAM: convert filtered BAM to SAM for downstream matrix reconstruction
 * 6) Convert_BAM_to_Counts: regenerate matrix.mtx using STARsolo features/barcodes plus filtered SAM
 * 7) Filter_raw_counts: apply STARsolo EmptyDrops_CR to drop empty barcodes and produce filtered matrices
 *
 * Inputs
 * ------
 * - FASTQs: `params.input` (expects *_R{1,2}_001.fastq.gz; paired-end)
 * - Genome resources: `params.fasta`, `params.fai`, `params.gtf`, `params.STAR_index`
 * - 10x barcode whitelist: `params.barcode_whitelist`
 * - WASP VCF for allele-aware alignment: `params.wasp_VCF`
 * - Helper scripts: `params.matrix_script1`, `params.matrix_script2`
 *
 * Outputs (key)
 * -------------
 * Per library (under `${params.outdir}/STAR/<sample_id>/`):
 *   - <sample_id>.Aligned.sortedByCoord.out.bam (+ index)
 *   - <sample_id>.Solo.out/ (STARsolo raw counts)
 *   - filtered_bams/ and reconstructed raw_counts/ + filtered_counts/
 *
 * Notes
 * -----
 * - Part 2 of this pipeline can re-use the STARsolo output directories produced here to perform
 *   demultiplexing (Souporcell) and genotype correlation.
 */
```


// -----------------------------
// Parameters
// -----------------------------
params.input             = '/gchm/scRNAseq/data/fastq_files/renamed/*_R{1,2}_001.fastq.gz'
params.outdir            = '/gchm/scRNAseq/sc_lib_preprocessing/results/03.27.24_preprocessing'
params.barcode_whitelist = '/gchm/scRNAseq/sc_lib_preprocessing/resources/10x_V3_barcode_whitelist.txt.gz'

params.aligner  = 'star'
params.protocol = '10XV3'
params.genome   = 'GRCh38'

// Variant resources
params.wasp_VCF    = '/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4p_WGS_pass_only.snps.maf1.WASP.vcf'
params.samples_VCF = '/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4_allsamples_common_maf1.vcf.gz'

// Genome resources (see 001_download_hg38.sh and 002_prepare_fasta.sh)
params.fasta      = '/gchm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa'
params.fai        = '/gchm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa.fai'
params.gtf        = '/gchm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf'
params.STAR_index = '/gchm/resources/genome/hg38_gencode_PRI_align'

// Helper scripts
params.matrix_script1 = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/soloCOUNTmatrixFromBam.py'
params.matrix_script2 = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/header_matrix.py'

params.Rscript   = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/batch_txt_script.r'
params.metadata  = '/gchm/scRNAseq/sc_lib_preprocessing/resources/cd4_scRNAseq_meta.csv'

params.expected_cells    = 25000
params.expected_clusters = 12

// -----------------------------
// Processes
// -----------------------------

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
    publishDir "${params.outdir}/STAR/${meta.sample_ID}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index_ch
    path wasp_VCF
    path barcode_whitelist
    val expected_cells

    output:
    tuple val(meta), path('*sortedByCoord.out.bam') , emit: bam
    tuple val(meta), path('*.Solo.out')             , emit: counts
    tuple val(meta), path('*Log.final.out')         , emit: log_final
    tuple val(meta), path('*Log.out')               , emit: log_out
    tuple val(meta), path('*Log.progress.out')      , emit: log_progress

    tuple val(meta), path('*toTranscriptome.out.bam'), optional: true, emit: bam_transcript
    tuple val(meta), path('*fastq.gz')               , optional: true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional: true, emit: tab

    script:
    // `reads` is a list of files for this sample; assumes R1/R2 present.
    def (forward, reverse) = reads.collate(2).transpose()

    """
    STAR --runMode alignReads \\
      --genomeDir ${index_ch} \\
      --readFilesIn ${reverse.join(",")} ${forward.join(",")} \\
      --outFileNamePrefix ${meta.sample_ID}. \\
      --soloCBwhitelist <(gzip -cdf ${barcode_whitelist}) \\
      --soloType CB_UMI_Simple \\
      --soloFeatures Gene \\
      --soloUMIlen 12 \\
      --outSAMtype BAM SortedByCoordinate \\
      --soloCellFilter EmptyDrops_CR ${expected_cells} 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
      --outSAMattrRGline ID:${meta.sample_ID} CN:NYGC SM:${meta.sample_ID} \\
      --readFilesCommand zcat \\
      --runDirPerm All_RWX \\
      --outWigType bedGraph \\
      --twopassMode Basic \\
      --waspOutputMode SAMtag \\
      --varVCFfile ${wasp_VCF} \\
      --limitOutSJcollapsed 5000000 \\
      --limitSjdbInsertNsj 1500000 \\
      --outSAMattributes NH HI AS NM MD vA vG vW GX CB UB
    """
}

process STAR_multiqc {
    publishDir "${params.outdir}/STAR_MultiQC/", mode: 'copy'

    input:
    tuple val(meta), path('*Log.final.out')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional: true, emit: plots

    script:
    """
    multiqc .
    """
}

process Bam_Filter {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
          path("*.sortedByCoord.filtered.out.bam"),
          path("*.sortedByCoord.filtered.out.bam.bai"),
          emit: bams

    script:
    """
    samtools index ${bam}

    # Keep reads that either lack vW tag or passed WASP (vW==1)
    samtools view -e '![vW] || [vW]==1' ${bam} -o ${meta.sample_ID}.temp.bam

    # Drop chrX/chrY (keeps header lines too)
    samtools view -h ${meta.sample_ID}.temp.bam \\
      | awk '{if(\$3 != "chrY" && \$3 != "chrX"){print \$0}}' \\
      | samtools view -Shb - > ${meta.sample_ID}.sortedByCoord.filtered.out.bam

    samtools index ${meta.sample_ID}.sortedByCoord.filtered.out.bam
    """
}

process BAM_to_SAM {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    """
    samtools view ${bam} > ${meta.sample_ID}.sortedByCoord.filtered.out.sam
    """
}

process Convert_BAM_to_Counts {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams/raw_counts", mode: 'copy'

    input:
    tuple val(meta), path(sam)
    tuple val(meta), path(counts)
    path matrix_script1
    path matrix_script2

    output:
    tuple val(meta), path("matrix.mtx")   , emit: matrix
    tuple val(meta), path("barcodes.tsv") , emit: barcodes
    tuple val(meta), path("features.tsv") , emit: features

    script:
    """
    cp ${counts}/Gene/raw/features.tsv features.tsv
    cp ${counts}/Gene/raw/barcodes.tsv barcodes.tsv

    python ${matrix_script1} ${counts}/Gene/raw/barcodes.tsv ${counts}/Gene/raw/features.tsv ${sam} \\
      | sort -k2,2n -k1,1n > matrix.mtx

    python ${matrix_script2} ${counts}/Gene/raw/features.tsv ${counts}/Gene/raw/barcodes.tsv matrix.mtx
    """
}

process Filter_raw_counts {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams/filtered_counts", mode: 'copy'

    input:
    tuple val(meta), path(matrix)
    tuple val(meta), path(barcodes)
    tuple val(meta), path(features)

    output:
    tuple val(meta), path("filtered.matrix.mtx")     , emit: matrix
    tuple val(meta), path("filtered.barcodes.tsv")   , emit: barcodes
    tuple val(meta), path("filtered.features.tsv")   , emit: features

    script:
    """
    STAR --runMode soloCellFiltering ./ filtered. \\
      --soloCellFilter EmptyDrops_CR 25000 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \\
      --runThreadN ${task.cpus}
    """
}

/*
 * Optional: Demultiplexing helpers (souporcell). Left here cleaned, but workflow calls are commented.
 */

process Demulti_prep_VCF {
    publishDir "${params.outdir}/resources", mode: 'copy'

    input:
    path samples_VCF

    output:
    path("CD4_allsamples_snps_maf1_sorted.vcf.gz"), emit: filtered_vcf
    path("*.bed")                                 , emit: exons_bed
    path("*.recode.vcf")                          , emit: exons_vars
    path("*.stats")                               , emit: exons_vars_stats

    script:
    """
    echo "initial stats"
    bcftools stats ${samples_VCF} > CD4_allsamples_common_maf1.stats

    echo "filtering for SNPs"
    bcftools view -v snps -Oz -o CD4_allsamples_snps_maf1.vcf.gz ${samples_VCF}
    bcftools sort CD4_allsamples_snps_maf1.vcf.gz -Oz -o CD4_allsamples_snps_maf1_sorted.vcf.gz
    bcftools stats CD4_allsamples_snps_maf1_sorted.vcf.gz > CD4_allsamples_snps_maf1_sorted.stats

    echo "dropping GT info"
    bcftools view -G -Oz -o CD4_snps_maf1.vcf.gz CD4_allsamples_snps_maf1_sorted.vcf.gz

    echo "downloading exon coordinates (hg38)"
    wget https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/_downloads/dbfbf640f118c0a8a98413c0731f9f9c/hg38exonsUCSC.bed

    awk '{print "chr"\$1, \$2, \$3, \$4, \$5, \$6}' hg38exonsUCSC.bed > formatted_hg38exonsUCSC.bed

    echo "filtering SNPs within exons"
    vcftools --gzvcf CD4_snps_maf1.vcf.gz \\
      --max-alleles 2 \\
      --remove-indels \\
      --bed formatted_hg38exonsUCSC.bed \\
      --recode \\
      --recode-INFO-all \\
      --out CD4_snps_maf1.exon_filtered

    bcftools stats CD4_snps_maf1.exon_filtered.recode.vcf > CD4_snps_maf1.exon_filtered.recode.vcf.stat
    """
}

process Mk_per_batch_txt {
    publishDir "${params.outdir}/Demultiplexing/perBatch_genotypes/sample_list", mode: 'copy'

    input:
    path Rscript
    path metadata

    output:
    path("*.txt"), emit: demulti_batch_txt

    script:
    """
    echo "Creating per-batch sample lists"
    Rscript ${Rscript} ${metadata}
    """
}

process Mk_per_batch_genotypes {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/Demultiplexing/perBatch_genotypes/batch_vcfs", mode: 'copy'

    input:
    tuple val(meta), path(txt_file)
    path filtered_vcf

    output:
    path("*.wgs_batch.vcf.gz"), emit: demulti_batch_vcfs

    script:
    """
    echo "Creating per-batch VCFs"
    bcftools view -S ${txt_file} --force-samples -Oz -o ${meta.sample_ID}.wgs_batch.vcf.gz ${filtered_vcf}
    """
}

process Demultiplexing_souporcell {
    tag "$meta.sample_ID"
    publishDir "${params.outdir}/Demultiplexing/souporcell", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(barcodes)
    val expected_clusters
    path fasta
    path exons_vars

    output:
    tuple val(meta), path("*.soupOut"), emit: demulti_batch

    script:
    """
    souporcell_pipeline.py \\
      -i ${bam} \\
      -b ${barcodes} \\
      -f ${fasta} \\
      -t ${task.cpus} \\
      -o "${meta.sample_ID}.soupOut" \\
      -k ${expected_clusters} \\
      --skip_remap SKIP_REMAP \\
      --common_variants ${exons_vars}
    """
}

// -----------------------------
// Workflow
// -----------------------------
workflow {

    Genome_Generate(params.fasta, params.fai, params.gtf, params.STAR_index)

    Channel
        .fromFilePairs(params.input)
        .map { id, reads ->
            def meta = [['sample_ID', 'library_id', 'lane', 'read', 'other'], id.tokenize("_")]
                .transpose()
                .collectEntries()
            [meta, reads]
        }
        // Keep only sample_ID in meta (so grouping works cleanly)
        .map { meta, reads -> [ meta.subMap('sample_ID'), reads ] }
        // Group multiple lanes per sample_ID
        .groupTuple(by: [0])
        .map { meta, reads -> [ meta, reads.flatten() ] }
        .set { group_fastqs }

    STAR_align(group_fastqs, Genome_Generate.out, params.wasp_VCF, params.barcode_whitelist, params.expected_cells)

    STAR_align.out.log_final
        | collectFile(name: '*Log.final.out')
        | set { multiqc_ch }

    // STAR_multiqc(multiqc_ch)

    Bam_Filter(STAR_align.out.bam)
    BAM_to_SAM(Bam_Filter.out.bams)

    Convert_BAM_to_Counts(
        BAM_to_SAM.out.sam,
        STAR_align.out.counts,
        params.matrix_script1,
        params.matrix_script2
    )

    Filter_raw_counts(
        Convert_BAM_to_Counts.out.matrix,
        Convert_BAM_to_Counts.out.barcodes,
        Convert_BAM_to_Counts.out.features
    )

    // Optional demux (kept here but not executed)
    // Demulti_prep_VCF(params.samples_VCF)
    // Mk_per_batch_txt(params.Rscript, params.metadata)
    // ...then Mk_per_batch_genotypes + Demultiplexing_souporcell
}
