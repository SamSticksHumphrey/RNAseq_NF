//
// Run RSeQC modules
//

include { RSEQC_BAMSTAT            } from '../modules/postAlignmentQC/rseqc/bamstat'
include { RSEQC_INNERDISTANCE      } from '../modules/postAlignmentQC/rseqc/innerdistance'
include { RSEQC_JUNCTIONANNOTATION } from '../modules/postAlignmentQC/rseqc/junctionannotation'
include { RSEQC_JUNCTIONSATURATION } from '../modules/postAlignmentQC/rseqc/junctionsaturation'
include { RSEQC_READDISTRIBUTION   } from '../modules/postAlignmentQC/rseqc/readdistribution'
include { RSEQC_READDUPLICATION    } from '../modules/postAlignmentQC/rseqc/readduplication'

workflow RSEQC {
    take:
    bam_bai       // channel: [ val(meta), [ bam_bai], [ bai ] ]
    bed           //    file: /path/to/genome.bed

    main:

    ch_versions = Channel.empty()

    //
    // Run RSeQC bam_stat.py
    //
    bamstat_txt = Channel.empty()

        RSEQC_BAMSTAT ( bam_bai)
        bamstat_txt = RSEQC_BAMSTAT.out.txt
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())


    //
    // Run RSeQC inner_distance.py
    //
    innerdistance_distance = Channel.empty()
    innerdistance_freq     = Channel.empty()
    innerdistance_mean     = Channel.empty()
    innerdistance_pdf      = Channel.empty()
    innerdistance_rscript  = Channel.empty()

        RSEQC_INNERDISTANCE ( bam_bai, bed )
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        ch_versions = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions.first())


    //
    // Run RSeQC junction_annotation.py
    //
    junctionannotation_bed          = Channel.empty()
    junctionannotation_interact_bed = Channel.empty()
    junctionannotation_xls          = Channel.empty()
    junctionannotation_pdf          = Channel.empty()
    junctionannotation_events_pdf   = Channel.empty()
    junctionannotation_rscript      = Channel.empty()
    junctionannotation_log          = Channel.empty()

        RSEQC_JUNCTIONANNOTATION ( bam_bai, bed )
        junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())


    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()
        RSEQC_JUNCTIONSATURATION ( bam_bai, bed )
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())


    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()
        RSEQC_READDISTRIBUTION ( bam_bai, bed )
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())


    //
    // Run RSeQC read_duplication.py
    //
    readduplication_seq_xls = Channel.empty()
    readduplication_pos_xls = Channel.empty()
    readduplication_pdf     = Channel.empty()
    readduplication_rscript = Channel.empty()

        RSEQC_READDUPLICATION ( bam_bai)
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())


    emit:
    bamstat_txt                     // channel: [ val(meta), txt ]

    innerdistance_distance          // channel: [ val(meta), txt ]
    innerdistance_freq              // channel: [ val(meta), txt ]
    innerdistance_mean              // channel: [ val(meta), txt ]
    innerdistance_pdf               // channel: [ val(meta), pdf ]
    innerdistance_rscript           // channel: [ val(meta), r   ]

    junctionannotation_bed          // channel: [ val(meta), bed ]
    junctionannotation_interact_bed // channel: [ val(meta), bed ]
    junctionannotation_xls          // channel: [ val(meta), xls ]
    junctionannotation_pdf          // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    junctionannotation_rscript      // channel: [ val(meta), r   ]
    junctionannotation_log          // channel: [ val(meta), log ]

    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    readduplication_seq_xls         // channel: [ val(meta), xls ]
    readduplication_pos_xls         // channel: [ val(meta), xls ]
    readduplication_pdf             // channel: [ val(meta), pdf ]
    readduplication_rscript         // channel: [ val(meta), r   ]

    versions = ch_versions          // channel: [ versions.yml ]
}
