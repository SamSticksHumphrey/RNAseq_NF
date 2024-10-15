process CAT_FASTQ {
    tag "$meta.sample_id"
    label 'process_low'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/fastq" },
        mode: 'copy',
        pattern: '*.fastq.gz',
        enabled: params.save_all_intermediates,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename },
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/fastq/checksum" },
        pattern: '*.fastq.gz.md5',
        mode: 'copy',
        enabled: params.save_all_intermediates
    ]]


    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    path "*.md5"                          , emit: md5
    path "versions.yml"                       , emit: versions

    script:
    def readList = reads.collect{ it.toString() }
    def read1 = []
    def read2 = []
    readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

    """
    cat ${read1.join(' ')} > ${meta.sample_id}_R1.fastq.gz
    cat ${read2.join(' ')} > ${meta.sample_id}_R2.fastq.gz

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch tmp.merged.fastq.gz
    touch versions.yml
    """
}
