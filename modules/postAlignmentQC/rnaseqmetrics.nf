//      description: |
//        A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS)
//        data and formats such as SAM/BAM/CRAM and VCF.
//      homepage: https://broadinstitute.github.io/picard/
//      documentation: https://broadinstitute.github.io/picard/
//      licence: ['MIT']
//

process PICARD_RNASEQMETRICS {
    tag "$meta.sample_id"
    label 'process_low'
    publishDir = [[
          path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/picard" },
          mode: 'copy',
          pattern: '*.rnaseq_metrics'
        ],[
          path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/picard/checksum" },
          mode: 'copy',
          pattern: '*.rnaseq_metrics.md5'
        ] ]


    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path refFlat
    path ribosomalIntervalList

    output:
    tuple val(meta), path("${meta.sample_id}.rnaseq_metrics") , emit: metrics
    path "*.md5"                                        , emit: md5
    path  "versions.yml"                                , emit: versions

    script:


    def avail_mem = 3
    if (!task.memory) {
        log.info '[CollectRnaSeqMetrics] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        CollectRnaSeqMetrics \\
        I=$bam \\
        O=${meta.sample_id}.rnaseq_metrics \\
        REF_FLAT=$refFlat \\
        STRAND=NONE \\
        RIBOSOMAL_INTERVALS=$ribosomalIntervalList

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
    touch ${meta.sample_id}.rnaseq_metrics
    touch versions.yml
    """
}
