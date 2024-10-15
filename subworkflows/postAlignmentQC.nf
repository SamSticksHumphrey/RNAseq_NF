//
// post alignment QC
//

include { SUBREAD_FEATURECOUNTS }       from '../modules/postAlignmentQC/subread_featurecounts'
include { MULTIQC_CUSTOM_BIOTYPE }      from '../modules/postAlignmentQC/multiqc_custom_biotype'
include { BEDTOOLS_GENOMECOV }          from '../modules/postAlignmentQC/bedtools_genomecov'
include { BEDCLIP          }            from '../modules/postAlignmentQC/bedclip'
include { BEDGRAPH_TO_BIGWIG }          from '../modules/postAlignmentQC/bedgraphtobigwig'
include { QUALIMAP_RNASEQ}              from '../modules/postAlignmentQC/qualimap'
include { DUPRADAR }                    from '../modules/postAlignmentQC/dupradar'
include { SAMTOOLS_MPILEUP }            from '../modules/postAlignmentQC/samtools_mpileup'
include { PICARD_RNASEQMETRICS }        from '../modules/postAlignmentQC/rnaseqmetrics'
include { RSEQC }                       from '../subworkflows/rseqc'

workflow POST_ALIGNMENT_QC {

  take:
    ch_genome_bam_bai
    ch_tx_bam
    ch_biotypes_header_multiqc
    fasta                       // channel: /path/to/fasta
    gtf
    ch_chrom_sizes
    ch_gene_bed
    ch_ref_flat
    ch_interval_list
    ch_ngscheckmate_snpbed

  main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run rnaseq metrics
    PICARD_RNASEQMETRICS ( ch_genome_bam_bai, ch_ref_flat, ch_interval_list )
    ch_versions = ch_versions.mix(PICARD_RNASEQMETRICS.out.versions.first())


    //
    // MODULE: Genome-wide coverage with rseqc
    RSEQC ( ch_genome_bam_bai, ch_gene_bed )
    ch_versions = ch_versions.mix(RSEQC.out.versions.first())


    ch_gtf = Channel.fromPath(gtf)

    //
    // MODULE: Feature biotype QC using featureCounts
    ch_gtf
      .map { biotypeInGtf(it, 'gene_biotype') }
      .set { biotype_in_gtf }

    // Prevent any samples from running if GTF file doesn't have a valid biotype
    ch_genome_bam_bai
        .combine ( ch_gtf )
        .combine ( biotype_in_gtf )
        .filter  { it[-1] }
        .map     { it[0..<it.size()-1] }
        .set     { ch_featurecounts }

    SUBREAD_FEATURECOUNTS ( ch_featurecounts )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    MULTIQC_CUSTOM_BIOTYPE ( SUBREAD_FEATURECOUNTS.out.counts, ch_biotypes_header_multiqc)
    ch_versions = ch_versions.mix(MULTIQC_CUSTOM_BIOTYPE.out.versions.first())

    //
    // MODULE: Genome-wide coverage with BEDTools
    BEDTOOLS_GENOMECOV (ch_genome_bam_bai)
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // Clip bedGraph file
    BEDCLIP ( BEDTOOLS_GENOMECOV.out.bedgraphs, ch_chrom_sizes )
    ch_versions = ch_versions.mix(BEDCLIP.out.versions.first())

    //
    // Convert bedGraph to bigWig
    BEDGRAPH_TO_BIGWIG ( BEDCLIP.out.bedClip_graphs, ch_chrom_sizes )
    ch_versions = ch_versions.mix(BEDGRAPH_TO_BIGWIG.out.versions.first())

    //
    // Run samtools mpileup over the NGSCheckMate SNP profile for sample
    //    provence checks later
    SAMTOOLS_MPILEUP( ch_genome_bam_bai, ch_ngscheckmate_snpbed, fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    //
    // MODULE: QUALIMAP
    QUALIMAP_RNASEQ ( ch_genome_bam_bai, gtf )
    ch_versions = ch_versions.mix(QUALIMAP_RNASEQ.out.versions.first())

    //
    // MODULE: DUPRADAR
    DUPRADAR ( ch_genome_bam_bai, gtf )
    ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())


    emit:

    picard_rnaseqmetrics = PICARD_RNASEQMETRICS.out.metrics

    dupradar_multiqc = DUPRADAR.out.multiqc


    featurecounts_multiqc = MULTIQC_CUSTOM_BIOTYPE.out.tsv

    bigwigs      = BEDGRAPH_TO_BIGWIG.out.bigwigs    // channel: [ val(meta), [ bigwig ] ]
    bedgraphs    = BEDCLIP.out.bedClip_graphs        // channel: [ val(meta), [ bedgraph ] ]

    qualimap_multiqc = QUALIMAP_RNASEQ.out.results


    ngscheckmate_vcf = SAMTOOLS_MPILEUP.out.vcffile

    bamstat_multiqc            = RSEQC.out.bamstat_txt
    innerdistance_multiqc      = RSEQC.out.innerdistance_freq
    junctionannotation_multiqc = RSEQC.out.junctionannotation_log
    junctionsaturation_multiqc = RSEQC.out.junctionsaturation_rscript
    readdistribution_multiqc   = RSEQC.out.readdistribution_txt
    readduplication_multiqc    = RSEQC.out.readduplication_pos_xls

    versions = ch_versions

}



//
// Function to check whether biotype field exists in GTF file
public static Boolean biotypeInGtf(gtf_file, biotype) {
   def hits = 0
   gtf_file.eachLine { line ->
       def attributes = line.split('\t')[-1].split()
       if (attributes.contains(biotype)) {
           hits += 1
       }
   }
   if (hits) {
       return true
   } else {
       log.info " Biotype attribute '${biotype}' not found"
       return false
   }
}
