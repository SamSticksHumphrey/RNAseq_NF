//
// Read QC, UMI extraction and trimming
//

include { OUTPUT_FASTQ     } from '../modules/preAlignmentQC/output_fastq'
include { CAT_FASTQ        } from '../modules/preAlignmentQC/cat_fastq'
include { FASTQC           } from '../modules/preAlignmentQC/fastqc'
include { TRIMGALORE       } from '../modules/preAlignmentQC/trimgalore'

workflow PRE_ALIGNMENT_QC {

  take:
  samplesheet   // file: /path/to/samplesheet.csv


  main:

    ch_versions = Channel.empty()

    ch_input = Channel.fromPath( samplesheet )

    //
    // MODULE: Check the samplesheet to ensure all files exist and formatted
    //  correctly. Calls check_samplesheet.py from bin/
    ch_input
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { ch_reads }

    //
    // MODULE: output all fastq files
    OUTPUT_FASTQ ( ch_reads )

    //
    // Take the output from SAMPLESHEET_CHECK and cat together any multi-fastq files
    ch_reads
        .groupTuple( by: [0] )
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ ( ch_fastq.multiple )
      .reads
      .mix(ch_fastq.single)
      .set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


    //
    // MODULE: run fastqc over all the the fastq files and output for multiqc
    //
    FASTQC ( ch_cat_fastq )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: run trimgalore over all the the fastq files and output
    //  trimming stats for multiqc
    //
    ch_trim_reads = ch_cat_fastq
    TRIMGALORE ( ch_cat_fastq ).reads.set { ch_trim_reads }
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())



    emit:
    reads = ch_trim_reads               // channel: [ val(meta), [ reads ] ]

    fastqc_html = FASTQC.out.html      // channel: [ val(meta), [ html ] ]
    fastqc_zip  = FASTQC.out.zip      // channel: [ val(meta), [ zip ] ]

    trim_html   = TRIMGALORE.out.html       // channel: [ val(meta), [ html ] ]
    trim_zip    = TRIMGALORE.out.zip       // channel: [ val(meta), [ zip ] ]
    trim_log    = TRIMGALORE.out.log      // channel: [ val(meta), [ txt ] ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}



// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.cellline_id         = row.cellline_id
    meta.sample_id           = row.sample_id

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    
    array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]

    return array
}
