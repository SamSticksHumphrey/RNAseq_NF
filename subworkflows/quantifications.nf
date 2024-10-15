//
// Gene/transcript quantifications
include { RSEM_CALCULATEEXPRESSION  } from '../modules/quantifications/rsem_calculateExpression'
include { KALLISTO_ALIGN            } from '../modules/quantifications/kallisto_align'
include { SALMON_QUANT              } from '../modules/quantifications/salmon_quantification'
include { SALMON_TXIMPORT           } from '../modules/quantifications/salmon_tximport'
include { RSEM_TXIMPORT             } from '../modules/quantifications/rsem_tximport'
include { FEATURECOUNTS_GENES_EXONS } from '../modules/quantifications/featurecounts_genes_exons'

workflow QUANTIFICATIONS {

    take:
    ch_reads
    ch_genome_bam_bai
    ch_tx_bam
    ch_transcript_fasta
    ch_gtf
    ch_bed
    ch_rsem_index
    ch_kallisto_index

    main:

    ch_versions = Channel.empty()

    //
    // Quantify reads with Kallisto
    KALLISTO_ALIGN ( ch_reads, ch_kallisto_index )
    ch_versions = ch_versions.mix(KALLISTO_ALIGN.out.versions.first())

    //
    // Quantify BAM with RSEM
    RSEM_CALCULATEEXPRESSION ( ch_tx_bam, ch_rsem_index )
    ch_versions = ch_versions.mix(RSEM_CALCULATEEXPRESSION.out.versions.first())

    RSEM_TXIMPORT ( RSEM_CALCULATEEXPRESSION.out.counts_transcript, ch_gtf )
    ch_versions = ch_versions.mix(RSEM_TXIMPORT.out.versions)

    //
    // Quantify and merge counts across samples
    SALMON_QUANT ( ch_tx_bam, ch_transcript_fasta, ch_gtf )
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())

    SALMON_TXIMPORT ( SALMON_QUANT.out.results, ch_gtf )
    ch_versions = ch_versions.mix(SALMON_TXIMPORT.out.versions)


    FEATURECOUNTS_GENES_EXONS ( ch_genome_bam_bai, ch_gtf )
    ch_versions = ch_versions.mix(FEATURECOUNTS_GENES_EXONS.out.versions.first())



    emit:

    rsem_stat                     = RSEM_CALCULATEEXPRESSION.out.stat              // channel: [ val(meta), stat ]
    rsem_counts_gene              = RSEM_CALCULATEEXPRESSION.out.counts_gene       // channel: [ val(meta), counts ]
    rsem_counts_transcript        = RSEM_CALCULATEEXPRESSION.out.counts_transcript // channel: [ val(meta), counts

    kallisto_abundace_h5          = KALLISTO_ALIGN.out.abundace_h5                 // channel: [ val(meta), HDF5  ]
    kallisto_abundace_tsv         = KALLISTO_ALIGN.out.abundace_tsv                // channel: [ val(meta), tsv   ]
    kallisto_run_info             = KALLISTO_ALIGN.out.run_info                    // channel: [ val(meta), json  ]
    kallisto_log                  = KALLISTO_ALIGN.out.log                    // channel: [ val(meta), json  ]

    salmon_results                       = SALMON_QUANT.out.results    // channel: [ val(meta), results_dir ]
    gene_tpm = SALMON_TXIMPORT.out.gene_tpms // path .tsv
    gene_counts = SALMON_TXIMPORT.out.gene_counts // path .tsv
    gene_lengths = SALMON_TXIMPORT.out.gene_lengths // path .tsv
    gene_counts_scaledByLibrarySizeAndGeneLength = SALMON_TXIMPORT.out.gene_counts_scaledByLibrarySizeAndGeneLength // path .tsv

    transcript_counts = SALMON_TXIMPORT.out.transcript_counts // path .tsv
    transcript_TPMs = SALMON_TXIMPORT.out.transcript_tpms // path .tsv
    transcript_lengths = SALMON_TXIMPORT.out.transcript_lengths // path .tsv
    transcript_scaledForDTU = SALMON_TXIMPORT.out.transcript_scaledForDTU // path .tsv

    featurecounts_multiqc = FEATURECOUNTS_GENES_EXONS.out.summary

    versions                 = ch_versions                                    // channel: [ versions.yml ]
}
