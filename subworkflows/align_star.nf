//
// Run start alignment unstranded so that rseqc can guess the strandedness
//

include { STAR_ALIGN        } from '../modules/align_star/star_align'
include { PICARD_MARKDUPLICATES } from '../modules/align_star/markduplicates'

include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_GENOME   } from '../modules/align_star/samtools_sort_index'
include { SAMTOOLS_SORT_INDEX as SAMTOOLS_SORT_INDEX_MARKDUP   } from '../modules/align_star/samtools_sort_index'
include { SAMTOOLS_ALLSTATS    } from '../modules/align_star/samtools_allstats'

workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    fasta // channel: /path/to/fasta
    gtf   // channel: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Align using STAR
    STAR_ALIGN ( reads, index, gtf )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    //
    // Sort Bam File into ??? format
    SAMTOOLS_SORT_INDEX_GENOME ( STAR_ALIGN.out.bam)

    //
    // Picard MarkDuplicates
    PICARD_MARKDUPLICATES (  SAMTOOLS_SORT_INDEX_GENOME.out.bam_bai )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    //
    // Sort Bam File into ??? format, and Index BAM file and run samtools stats, flagstat and idxstats
    SAMTOOLS_SORT_INDEX_MARKDUP ( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_MARKDUP.out.versions.first())

    //
    // Run samtools stats, flagstat and idxstats
    SAMTOOLS_ALLSTATS (SAMTOOLS_SORT_INDEX_MARKDUP.out.bam_bai, fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_ALLSTATS.out.versions.first())


    emit:
    orig_bam         = STAR_ALIGN.out.bam             // channel: [ val(meta), bam            ]
    log_final        = STAR_ALIGN.out.log_final       // channel: [ val(meta), log_final      ]
    log_out          = STAR_ALIGN.out.log_out         // channel: [ val(meta), log_out        ]
    log_progress     = STAR_ALIGN.out.log_progress    // channel: [ val(meta), log_progress   ]
    tx_bam           = STAR_ALIGN.out.bam_transcript  // channel: [ val(meta), bam_transcript ]

    tab              = STAR_ALIGN.out.tab             // channel: [ val(meta), tab            ]

    metrics          = PICARD_MARKDUPLICATES.out.metrics   // channel: [ val(meta), [ metrics ] ]
    genome_bam_bai   = SAMTOOLS_SORT_INDEX_MARKDUP.out.bam_bai        // channel: [ val(meta), [ bam ], [bai]]
    stats            = SAMTOOLS_ALLSTATS.out.stats            // channel: [ val(meta), [ stats ] ]
    flagstat         = SAMTOOLS_ALLSTATS.out.flagstat      // channel: [ val(meta), [ flagstat ] ]
    idxstats         = SAMTOOLS_ALLSTATS.out.idxstats      // channel: [ val(meta), [ idxstats ] ]

    versions         = ch_versions                    // channel: [ versions.yml ]
}
