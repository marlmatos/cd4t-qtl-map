#!/usr/bin/env nextflow

/*
 * scRNA-seq preprocessing pipeline (Part 2): post-alignment filtering, matrix rebuild, and genotype demultiplexing
 *
 * Based on: nf-core/scrnaseq (adapted)
 * Modified by: Marliette Matos
 * Date: <2024-05-24>
 *
 * Purpose
 * -------
 * This workflow continues from Part 1 outputs. It assumes STARsolo alignment results already exist on disk
 * (i.e., per-library STAR output directories containing Aligned.sortedByCoord.out.bam and Solo.out/).
 *
 * Major steps
 * -----------
 * 1) Load existing STARsolo outputs from `params.STAR_align_output`
 * 2) Filter aligned BAMs to retain autosomal reads and WASP-passing reads (vW==1 or missing vW tag)
 * 3) Rebuild raw Gene expression count matrix (matrix.mtx) from the filtered BAM/SAM plus STARsolo features/barcodes
 * 4) Apply STARsolo EmptyDrops-based cell filtering to generate filtered matrices
 * 5) Prepare a common-variant VCF for demultiplexing (exonic SNP subset)
 * 6) Run Souporcell demultiplexing per library
 * 7) Generate per-batch WGS VCFs and correlate Souporcell clusters to donor genotypes
 *
 * Inputs expected from Part 1
 * ---------------------------
 * For each library directory matched by `params.STAR_align_output` (e.g., .../STAR/<SAMPLE_ID>/):
 *   - <SAMPLE_ID>.Aligned.sortedByCoord.out.bam
 *   - <SAMPLE_ID>.Solo.out/Gene/raw/barcodes.tsv
 *   - <SAMPLE_ID>.Solo.out/Gene/raw/features.tsv
 *
 * Notes
 * -----
 * - This script does NOT run STAR alignment; alignment is performed in Part 1.
 * - Dependencies are filesystem-based (pt2 reads pt1 outputs from disk), not Nextflow channel-based.

* Parts of this pipeline (cluster-to-donor assignment via genotype correlation and related helper scripts)
 * are adapted from Demuxafy. We modified the original code for our study-specific inputs/outputs and
 * integration into this Nextflow workflow.
 */
// -----------------------------------------------------------------------------
// Parameters / directories
// -----------------------------------------------------------------------------

params.input              = '/gchm/scRNAseq/data/fastq_files/renamed/*_R{1,2}_001.fastq.gz'
params.outdir             = '/gchm/scRNAseq/sc_lib_preprocessing/results/03.27.24_preprocessing'
params.barcode_whitelist  = '/gchm/scRNAseq/sc_lib_preprocessing/resources/10x_V3_barcode_whitelist.txt.gz'

params.aligner            = 'star'
params.protocol           = '10XV3'
params.genome             = 'GRCh38'

params.wasp_VCF           = '/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4p_WGS_pass_only.snps.maf1.WASP.vcf'      // created with 001_preparing_wasp_vcf.sh
params.samples_VCF        = '/gchm/scRNAseq/sc_lib_preprocessing/resources/variant_vcf/CD4_allsamples_common_maf1.vcf.gz'            // created with 002_preparing_demultiplexing_vcf.sh
params.commmon_vars       = '/gchm/reference_genomes/1KG_common_SNPs/common_variants_grch38_chr_format.vcf'

// see 001_download_hg38.sh and 002_prepare_fasta.sh
params.fasta              = '/gchm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa'
params.fai                = '/gchm/resources/genome/hg38_gencode_PRI_align/GRCh38.primary_assembly_subset_masked.genome.fa.fai'
params.gtf                = '/gchm/resources/genome/hg38_gencode_raw/gencode.v44.primary_assembly.annotation.gtf'
params.STAR_index         = '/gchm/resources/genome/hg38_gencode_PRI_align'

params.matrix_script1     = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/soloCOUNTmatrixFromBam.py'
params.matrix_script2     = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/header_matrix.py'
params.corr_script        = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/demulti_corr.r'

params.Rscript            = '/gchm/scRNAseq/sc_lib_preprocessing/scripts/scrnaseq_preprocessing_nextflow/batch_txt_script.r'
params.metadata           = '/gchm/metadata_harmonizing/results/post_WGS_QC/version_080923/revision_103023/cd4_scRNAseq_meta.v2.csv'

params.STAR_align_output  = '/gchm/scRNAseq/sc_lib_preprocessing/results/03.27.24_preprocessing/STAR/scGEX-A*'

params.expected_cells     = 25000
params.expected_clusters  = 12


// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------

// Step 1: Create STAR index file
process Genome_Generate {

    publishDir "${params.STAR_index}", mode: 'copy'

    input:
    path(fasta)
    path(fai)
    path(gtf)
    path(STAR_index)

    output:
    path('STAR')

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


/*
 * Step 2: Align scRNAseq reads to the reference genome hg38
 *
 * Using STAR-solo, a VCF with genotype variants with MAF >0.01 is provided
 * for allele specific alignment using --WASPOutMode.
 *
 * SAMTags "GX CB UB" need to be added to the BAM file to recreate the matrix.mtx
 */
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
    def (forward, reverse) = reads.collate(2).transpose()

    """
    STAR --runMode alignReads \\
      --genomeDir ${index_ch} \\
      --readFilesIn ${reverse.join(",")} ${forward.join(",")} \\
      --outFileNamePrefix ${meta.sample_ID}. \\
      --soloCBwhitelist <(gzip -cdf $barcode_whitelist) \\
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


/*
 * Step 2.1: Work-around to provide the pipeline the output of STAR_align
 *
 * The job was killed before the last two libraries finished, so the output
 * stayed in the temp directory. In a time crunch, the completed output from
 * the libraries that were missing A9 and A29 was copied to the publish dir.
 *
 * A "touch process" is used to provide the outputs via a channel.
 */


/*
 * Step 3: Post-alignment QC summary
 */
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


/*
 * Step 4: Aligned BAM filtering
 *
 * Since we are using the WASPoutMode from STAR, STAR adds the flag "vW" to
 * the BAM file, but it does not take it into account for filtering or creating
 * the matrix.mtx.
 *
 * Here we filter for reads that PASSED WASP filtering (vW==1) and reads that do
 * not contain the WASP tag because no reads colocalize with the variant.
 *
 * Also filtering for reads mapping only to autosomes for downstream analysis.
 */
process Bam_Filter {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(paths)

    output:
    tuple val(meta), path("*.sortedByCoord.filtered.out.bam"), path("*.sortedByCoord.filtered.out.bam.bai"), emit: bams

    script:
    """
    samtools index ${paths}/${meta.sample_ID}.Aligned.sortedByCoord.out.bam

    samtools view -e '![vW] || [vW]==1' \\
      ${paths}/${meta.sample_ID}.Aligned.sortedByCoord.out.bam \\
      -o ${meta.sample_ID}.temp.bam

    samtools view -h ${meta.sample_ID}.temp.bam \\
      | awk '{if(\$3 != "chrY" && \$3 != "chrX"){print \$0}}' \\
      | samtools view -Shb - \\
      > ${meta.sample_ID}.sortedByCoord.filtered.out.bam

    samtools index ${meta.sample_ID}.sortedByCoord.filtered.out.bam
    """
}


/*
 * Step 5: Convert BAM to SAM needed for downstream analysis
 */
process BAM_to_SAM {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    """
    samtools view $bam > ${meta.sample_ID}.sortedByCoord.filtered.out.sam
    """
}


/*
 * Step 6: Recalculate 'matrix.mtx' from the filtered bam file
 */
process Convert_BAM_to_Counts {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams/raw_counts", mode: 'copy'

    input:
    tuple val(meta), path(sam), path(paths)
    path(matrix_script1)
    path(matrix_script2)

    output:
    tuple val(meta), path("matrix.mtx"), path("barcodes.tsv"), path("features.tsv")

    script:
    """
    cp $paths/${meta.sample_ID}.Solo.out/Gene/raw/features.tsv features.tsv
    cp $paths/${meta.sample_ID}.Solo.out/Gene/raw/barcodes.tsv barcodes.tsv

    python $matrix_script1 \\
      $paths/${meta.sample_ID}.Solo.out/Gene/raw/barcodes.tsv \\
      $paths/${meta.sample_ID}.Solo.out/Gene/raw/features.tsv \\
      $sam \\
      | sort -k2,2n -k1,1n \\
      > matrix.mtx

    python $matrix_script2 \\
      $paths/${meta.sample_ID}.Solo.out/Gene/raw/features.tsv \\
      $paths/${meta.sample_ID}.Solo.out/Gene/raw/barcodes.tsv \\
      matrix.mtx
    """
}


/*
 * Step 6: Filter out empty barcodes (99 percentile) from the recalculated bam files
 */
process Filter_raw_counts {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/STAR/${meta.sample_ID}/filtered_bams/filtered_counts", mode: 'copy'

    input:
    tuple val(meta), path(matrix), path(barcodes), path(features)

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
 * Step 7: Genotype demultiplexing of scRNAseq batches
 *
 * In the experimental design for this scRNAseq experiment, CD4 T cells from
 * 408 distinct donors were pooled in batches (12 donors each), for a total
 * of 34 single cell 10x GEX libraries.
 *
 * The reads pertaining each sample can be demultiplexed from each batch by
 * associating short read variant calling (from scRNAseq) with WGS sequencing
 * of each sample.
 *
 * Software: Souporcell
 * https://www.nature.com/articles/s41592-020-0820-1
 */


/*
 * Step 7.1: Prepare the reference VCF with the variants needed to demultiplex scRNA-seq experiments.
 *
 * - Gather the reference SNP VCF file and filter for MAF>1 or MAF>5
 * - Drop genotype "GT" for all samples, only need variant position
 * - Downloading exonic coordinates for hg38
 * - The hg38exonsUCSC.bed file does not have the same "chr" formatting as the CD4 SNPs VCF file; add it
 * - Only keep SNPs within exonic regions
 *
 * Note: using the "Demuxafy.sif" singularity image because it already contains all packages needed
 */
process Demulti_prep_VCF {

    publishDir "${params.outdir}/resources", mode: 'copy'

    input:
    path(samples_VCF)

    output:
    path("CD4_allsamples_snps_maf1.vcf.gz") , emit: filtered_vcf
    path("*.bed")                          , emit: exons_bed
    path("*.recode.vcf")                   , emit: exons_vars
    path("*.stats")                        , emit: exons_vars_stats

    script:
    """
    echo "initial stats"
    bcftools stats $samples_VCF > CD4_allsamples_common_maf1.stats

    echo "filtering for SNPS"
    bcftools view -v snps -Oz -o CD4_allsamples_snps_maf1.vcf.gz $samples_VCF

    bcftools stats CD4_allsamples_snps_maf1.vcf.gz > CD4_allsamples_snps_maf1.stats

    echo "Dropping genotype (GT) info for all samples"
    bcftools view -G -Oz -o CD4_snps_maf1.vcf.gz CD4_allsamples_snps_maf1.vcf.gz

    echo "downloading exonic coordinates for hg38"
    wget https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/_downloads/dbfbf640f118c0a8a98413c0731f9f9c/hg38exonsUCSC.bed

    awk '{print "chr"\$1, \$2, \$3, \$4, \$5, \$6}' hg38exonsUCSC.bed > formatted_hg38exonsUCSC.bed

    echo "filtering for SNPs within exonic regions only"
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


/*
 * Step 7.2: Run demultiplexing of scRNAseq batches
 *
 * Souporcell SINGULARITY container is at:
 * https://singularityhub.github.io/singularityhub-archive/containers/wheaton5-souporcell-latest/
 *
 * It was downloaded to the local machine and with nextflow, see nextflow.config file.
 */
process Demultiplexing_souporcell {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/Demultiplexing/souporcell", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes)
    val(expected_clusters)
    path(fasta)
    path(exons_vars)

    output:
    tuple val(meta), path("*.soupOut"), emit: demulti_batch

    script:
    """
    souporcell_pipeline.py \\
      -i $bam \\
      -b $barcodes \\
      -f $fasta \\
      -t ${task.cpus} \\
      -o "${meta.sample_ID}.soupOut" \\
      -k $expected_clusters \\
      --skip_remap SKIP_REMAP \\
      --common_variants $exons_vars
    """
}


/*
 * Step 7.3: Prepare .txt files with the sample names in each batch as written in the reference SNP file
 *
 * Rationale:
 * Since souporcell demultiplexing uses common variants, correlate the real
 * genotypes from each batch of WGS with the genotypes from the scRNAseq experiment.
 *
 * Creates txt files with sample names in each genotyping batch as written in
 * the metadata. Used to create batch-specific VCF files with the 12 samples in
 * each scRNAseq batch. Correlate with Pearson correlation.
 *
 * Note: some samples were not genotyped and appear as NA in the batches.
 */
process Mk_per_batch_txt {

    publishDir "${params.outdir}/Demultiplexing/perBatch_genotypes/sample_list", mode: 'copy'

    input:
    path(Rscript)
    path(metadata)

    output:
    path("*.txt"), emit: demulti_batch_txt

    script:
    """
    echo "Creating .txt files with the sample names in each genotyping batch"
    Rscript $Rscript $metadata
    """
}


/*
 * Step 7.4: Split the VCF containing the genotypes of all samples into 34 VCFs
 * containing only the samples of each scRNAseq batch.
 */
process Mk_per_batch_genotypes {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/Demultiplexing/perBatch_genotypes/batch_vcfs", mode: 'copy'

    input:
    tuple val(meta), path(txt_file)
    path(filtered_vcf)

    output:
    tuple val(meta), path("*.wgs_batch.vcf.gz"), emit: demulti_batch_vcfs

    script:
    """
    echo "Creating .vcf files for each genotyping batch"
    bcftools view -S $txt_file --force-samples -Oz -o ${meta.sample_ID}.wgs_batch.vcf.gz $filtered_vcf
    """
}


/*
 * Step 7.5: Correlate the demultiplexed clusters from Souporcell with the real
 * genotypes of each sample.
 */
process Cluster_genotype_corr {

    tag "$meta.sample_ID"
    publishDir "${params.outdir}/Demultiplexing//souporcell/${meta.sample_ID}.soupOut", mode: 'copy'

    input:
    tuple val(meta), path(cluster_genotypes), path(Ref_VCF)
    path(corr_script)

    output:
    tuple val(meta),
          path("*_key.txt"),
          path("*_correlation.png"),
          path("*_correlations.tsv"),
          emit: corr_batch

    script:
    """
    echo "Correlacting demutiplexed clusters from '$meta.sample_ID' Souporcell with a the real genotypes of each sample"
    cp $cluster_genotypes/cluster_genotypes.vcf cluster_genotypes.vcf

    $corr_script -r $Ref_VCF -c cluster_genotypes.vcf -o .
    """
}


// Rscript $corr_script $Ref_VCF $cluster_genotypes/common_variants_covered.vcf "${meta.sample_ID}.corrOut"


// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------

workflow {

    // Genome_Generate(params.fasta, params.fai, params.gtf, params.STAR_index)

    // Channel.fromFilePairs(params.input)
    //   | map { id, reads ->
    //       meta = [['sample_ID', 'library_id', 'lane', 'read', 'other'], id.tokenize("_")]
    //           .transpose()
    //           .collectEntries()
    //       [meta, reads]
    //   }
    //   | map { meta, reads -> [meta.subMap('sample_ID'), reads] }
    //   | groupTuple(by: [0])
    //   | map { meta, reads -> [meta, reads.flatten()] }
    //   | set { group_fastqs }

    // STAR_align(group_fastqs, Genome_Generate.out, params.wasp_VCF, params.barcode_whitelist, params.expected_cells)

    // STAR_align.out.log_final
    //   | collectFile(name: '*Log.final.out')
    //   | set { multiqc_ch }

    Channel.fromPath(params.STAR_align_output, type: 'dir')
        | collect
        | map { path -> [path.baseName - path.extension, path] }
        | transpose
        | map { id, path ->
            meta = [['sample_ID'], id.tokenize("")]
                .transpose()
                .collectEntries()
            [meta, path]
        }
        | view()
        | set { touch_align_ch }

    // STAR_multiqc(multiqc_ch)

    Bam_Filter(touch_align_ch)
    BAM_to_SAM(Bam_Filter.out.bams)

    // Combine multiple output channels so the libraries coincide, joining them by meta:sample_ID
    BAM_to_SAM.out.sam.combine(touch_align_ch, by: 0).view().set { convert_to_counts_input_ch }

    Convert_BAM_to_Counts(convert_to_counts_input_ch, params.matrix_script1, params.matrix_script2)

    // Pass the combined channel to Filter_raw_counts
    Filter_raw_counts(Convert_BAM_to_Counts.out)

    // Combine multiple output channels so the libraries coincide, joining them by meta:sample_ID
    Bam_Filter.out.bams.combine(Filter_raw_counts.out.barcodes, by: 0).view().set { demulti_input }

    Demulti_prep_VCF(params.samples_VCF)

    Demultiplexing_souporcell(
        demulti_input,
        params.expected_clusters,
        params.fasta,
        Demulti_prep_VCF.out.exons_vars
    )

    Mk_per_batch_txt(params.Rscript, params.metadata)

    Mk_per_batch_txt.out.demulti_batch_txt
        | collect
        | map { path -> [path.baseName - path.extension, path] }
        | transpose
        | map { id, path ->
            meta = [['sample_ID'], id.tokenize("")]
                .transpose()
                .collectEntries()
            [meta, path]
        }
        | view
        | set { wgs_batches_ch }

    Mk_per_batch_genotypes(wgs_batches_ch, Demulti_prep_VCF.out.filtered_vcf)

    Demultiplexing_souporcell.out.demulti_batch
        .combine(Mk_per_batch_genotypes.out.demulti_batch_vcfs, by: 0)
        .view()
        .set { corr_input_ch }

    Cluster_genotype_corr(corr_input_ch, params.corr_script)
}
