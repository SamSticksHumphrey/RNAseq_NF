#!/usr/bin/env nextflow

/*
========================================================================================
RNA-seq pipeline
========================================================================================
*/

nextflow.enable.dsl = 2

// Check input path parameters to see if they exist
checkPathParamList = [
      params.input, params.genome_fasta, params.transcript_fasta, params.genome_gtf,
      params.star_index, params.kallisto_index, params.rsem_index, params.genome_bed, params.genome_sizes,
      params.ngscheckmate_snpbed, params.ribosomal_interval_list, params.ref_flat
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Add a batchDir output

/*
========================================================================================
    Config and reference files
========================================================================================
*/

ch_multiqc_config            = file("$projectDir/configurationAndFormatting/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_individual_config = file("$projectDir/configurationAndFormatting/individual_multiqc_config.yaml", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/configurationAndFormatting/biotypes_header.txt", checkIfExists: true)

ch_fasta = file(params.genome_fasta)
ch_gtf = file(params.genome_gtf)
ch_bed = file(params.genome_bed)
ch_sizes = file(params.genome_sizes)
ch_transcript_fasta = file(params.transcript_fasta)
ch_star_index = file(params.star_index)
ch_rsem_index = file(params.rsem_index)
ch_kallisto_index = file(params.kallisto_index)
ch_ref_flat   = file(params.ref_flat)
ch_interval_list   = file(params.ribosomal_interval_list)
ch_ngscheckmate_snpbed = file(params.ngscheckmate_snpbed)


/*
================================================================================
    Import modules and subworkflows
================================================================================
*/

// PRE_ALIGNMENT_QC
include{ PRE_ALIGNMENT_QC } from "./subworkflows/preAlignmentQC"

// ALIGNMENTS
include{ ALIGN_STAR } from "./subworkflows/align_star"
include{ QUANTIFICATIONS } from "./subworkflows/quantifications"

// POST_ALIGNMENT_QC
include{ POST_ALIGNMENT_QC } from "./subworkflows/postAlignmentQC"

// MODULES
include { RSEQC_INFEREXPERIMENT }       from './modules/postAlignmentQC/rseqc/inferexperiment'
include { RSEQC_TIN             }       from './modules/postAlignmentQC/rseqc/tin'
include { RUN_NGSCHECKMATE      }       from "./modules/ngscheckmate_run"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "./modules/dumpsoftwareversions"
include { MULTIQC_BATCH         }       from "./modules/multiqc_batch"
include { MULTIQC_INDIVIDUAL    }       from "./modules/multiqc_individual"


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, run qc checks and trimming
    //  Currently ribosomal RNA is not removed from main bam file, but output as
    //    visualisation
     PRE_ALIGNMENT_QC (
        ch_input
    )
    ch_reads_fastq = PRE_ALIGNMENT_QC.out.reads
    ch_versions = ch_versions.mix( PRE_ALIGNMENT_QC.out.versions)


    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    ALIGN_STAR (
        ch_reads_fastq,
        ch_star_index,
        ch_fasta,
        ch_gtf
    )
    ch_genome_bam_bai   = ALIGN_STAR.out.genome_bam_bai
    ch_tx_bam           = ALIGN_STAR.out.tx_bam
    ch_versions         = ch_versions.mix(ALIGN_STAR.out.versions)


    //
    // Run RSeQC tin.py
    // This process takes a long time to complete since it only runs on 1 processor. Included here to get started ASAP
    tin_multiqc = ch_tx_bam.map{meta, bam -> [meta, []]}

    if(!params.skip_rseqc_tin) {
        RSEQC_TIN ( ch_genome_bam_bai, ch_bed)
        tin_multiqc     = RSEQC_TIN.out.tin_multiqc
        ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())
      }


    //
    // Inferexperiment for strandedness
    RSEQC_INFEREXPERIMENT( ch_genome_bam_bai, ch_bed )
    ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions)

    RSEQC_INFEREXPERIMENT.out.inferexperiment_txt
            .map { meta, strand_log -> [meta, getInferexperimentStrandedness(strand_log, 30)] }
            .set { ch_strand }


    // Add strand inormation to the bam and fastq files
    ch_reads_fastq.join(ch_strand, by: 0).set{ch_reads_fastq}
    ch_genome_bam_bai.join(ch_strand, by: 0).set{ch_genome_bam_bai}
    ch_tx_bam.join(ch_strand, by: 0).set{ch_tx_bam}

    //
    // SUBWORKFLOW: Count reads from BAM alignments using Salmon
    QUANTIFICATIONS (
      ch_reads_fastq,
      ch_genome_bam_bai,
      ch_tx_bam,
      ch_transcript_fasta,
      ch_gtf,
      ch_bed,
      ch_rsem_index,
      ch_kallisto_index
    )
    ch_versions = ch_versions.mix(QUANTIFICATIONS.out.versions)


    //
    // SUBWORKFLOW: post-alignment QC checks
    POST_ALIGNMENT_QC(
      ch_genome_bam_bai,
      ch_tx_bam,
      ch_biotypes_header_multiqc,
      ch_fasta,
      ch_gtf,
      ch_sizes,
      ch_bed,
      ch_ref_flat,
      ch_interval_list,
      ch_ngscheckmate_snpbed
    )
    ch_versions = ch_versions.mix(POST_ALIGNMENT_QC.out.versions)

    ch_ngscheckmate_vcfs = POST_ALIGNMENT_QC.out.ngscheckmate_vcf

    RUN_NGSCHECKMATE(ch_ngscheckmate_vcfs.collect(), ch_ngscheckmate_snpbed)

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_versions_mqc_yml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.first()


    // To gather all QC reports for MultiQC
    ch_report_individual = ALIGN_STAR.out.log_final
                            .join(PRE_ALIGNMENT_QC.out.fastqc_zip, by: 0)
                            .join(PRE_ALIGNMENT_QC.out.trim_zip, by: 0)
                            .join(PRE_ALIGNMENT_QC.out.trim_log, by: 0)
                            .join(ALIGN_STAR.out.metrics, by: 0)
                            .join(POST_ALIGNMENT_QC.out.picard_rnaseqmetrics, by: 0)
                            .join(ALIGN_STAR.out.stats, by: 0)
                            .join(ALIGN_STAR.out.flagstat, by: 0)
                            .join(ALIGN_STAR.out.idxstats, by: 0)
                            .join(QUANTIFICATIONS.out.rsem_stat, by: 0)
                            .join(QUANTIFICATIONS.out.salmon_results, by: 0)
                            .join(QUANTIFICATIONS.out.kallisto_log, by: 0)
                            .join(QUANTIFICATIONS.out.featurecounts_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.dupradar_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.featurecounts_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.qualimap_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.bamstat_multiqc, by: 0)
                            .join(RSEQC_INFEREXPERIMENT.out.inferexperiment_txt, by: 0)
                            .join(POST_ALIGNMENT_QC.out.innerdistance_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.junctionannotation_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.junctionsaturation_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.readdistribution_multiqc, by: 0)
                            .join(POST_ALIGNMENT_QC.out.readduplication_multiqc, by: 0)
                            .join(tin_multiqc, by: 0)

    //
    // MODULE: MultiQC
    MULTIQC_INDIVIDUAL (
        ch_multiqc_individual_config,
        ch_versions_mqc_yml.ifEmpty([]),
        ch_report_individual
        )

    multiqc_ind_report = MULTIQC_INDIVIDUAL.out.report.toList()

    // To gather all QC reports for MultiQC
    ch_report_batch  = Channel.empty()
            .mix(PRE_ALIGNMENT_QC.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
            .mix(PRE_ALIGNMENT_QC.out.trim_zip.collect{it[1]}.ifEmpty([]))
            .mix(PRE_ALIGNMENT_QC.out.trim_log.collect{it[1]}.ifEmpty([]))
            .mix(ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
            .mix(ALIGN_STAR.out.metrics.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.picard_rnaseqmetrics.collect{it[1]}.ifEmpty([]))
            .mix(ALIGN_STAR.out.stats.collect{it[1]}.ifEmpty([]))
            .mix(ALIGN_STAR.out.flagstat.collect{it[1]}.ifEmpty([]))
            .mix(ALIGN_STAR.out.idxstats.collect{it[1]}.ifEmpty([]))
            .mix(QUANTIFICATIONS.out.rsem_stat.collect{it[1]}.ifEmpty([]))
            .mix(QUANTIFICATIONS.out.salmon_results.collect{it[1]}.ifEmpty([]))
            .mix(QUANTIFICATIONS.out.kallisto_log.collect{it[1]}.ifEmpty([]))
            .mix(QUANTIFICATIONS.out.featurecounts_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.dupradar_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.featurecounts_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.qualimap_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.bamstat_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(RSEQC_INFEREXPERIMENT.out.inferexperiment_txt.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.innerdistance_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.junctionannotation_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.readdistribution_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(POST_ALIGNMENT_QC.out.readduplication_multiqc.collect{it[1]}.ifEmpty([]))
            .mix(tin_multiqc.collect{it[1]}.ifEmpty([]))


    //
    // MODULE: MultiQC
    MULTIQC_BATCH (
        ch_multiqc_config,
        ch_versions_mqc_yml.ifEmpty([]),
        ch_report_batch.collect()
        )
    multiqc_report = MULTIQC_BATCH.out.report.toList()

}

//
// Function that parses and returns the alignment rate from the STAR log output
//
public static ArrayList getStarPercentMapped(params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    def pass = false
    if (percent_aligned >= params.min_mapped_reads.toFloat()) {
        pass = true
    }
    return [ percent_aligned, pass ]
}


//
// Function that parses and returns the predicted strandedness from the RSeQC infer_experiment.py output
//
static String getInferexperimentStrandedness(inferexperiment_file, cutoff=30) {
    def sense        = 0
    def antisense    = 0
    def undetermined = 0
    inferexperiment_file.eachLine { line ->
        def undetermined_matcher = line =~ /Fraction of reads failed to determine:\s([\d\.]+)/
        def se_sense_matcher     = line =~ /Fraction of reads explained by "\++,--":\s([\d\.]+)/
        def se_antisense_matcher = line =~ /Fraction of reads explained by "\+-,-\+":\s([\d\.]+)/
        def pe_sense_matcher     = line =~ /Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s([\d\.]+)/
        def pe_antisense_matcher = line =~ /Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s([\d\.]+)/
        if (undetermined_matcher) undetermined = undetermined_matcher[0][1].toFloat() * 100
        if (se_sense_matcher)     sense        = se_sense_matcher[0][1].toFloat() * 100
        if (se_antisense_matcher) antisense    = se_antisense_matcher[0][1].toFloat() * 100
        if (pe_sense_matcher)     sense        = pe_sense_matcher[0][1].toFloat() * 100
        if (pe_antisense_matcher) antisense    = pe_antisense_matcher[0][1].toFloat() * 100
    }
    def strandedness = 'unstranded'
    if (sense >= 100-cutoff) {
        strandedness = 'forward'
    } else if (antisense >= 100-cutoff) {
        strandedness = 'reverse'
    }
    return strandedness
}

/*
========================================================================================
    THE END
========================================================================================
*/
