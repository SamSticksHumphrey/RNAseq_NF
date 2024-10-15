//      description: |
//        A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS)
//        data and formats such as SAM/BAM/CRAM and VCF.
//      homepage: https://broadinstitute.github.io/picard/
//      documentation: https://broadinstitute.github.io/picard/
//      licence: ['MIT']
//

process PICARD_MARKDUPLICATES {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/picard" },
            mode: 'copy',
            pattern: '*metrics.txt'
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/picard/checksum" },
            mode: 'copy',
            pattern: '*metrics.txt.md5'
        ]]


    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam")                  , emit: bam
    tuple val(meta), path("*.metrics.txt")          , emit: metrics
    path  "versions.yml"                            , emit: versions
    path "*.md5"                                    , emit: md5

    script:
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga.intdiv(2)
    }

    """
    picard \\
        -Xmx${avail_mem}g \\
        MarkDuplicates \\
        ASSUME_SORTED=true \\
        REMOVE_DUPLICATES=false \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp \\
        I=$bam \\
        O=${meta.sample_id}_Aligned_markdup.bam \\
        M=${meta.sample_id}_Aligned_markdup.metrics.txt


    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    """
    touch tmp_markduplicates.bam
    touch tmp.metrics.txt
    touch versions.yml
    """
}
