process OUTPUT_FASTQ {
    tag "$meta.sample_id"
    stageInMode 'copy'
    label 'process_low'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/fastq" },
        mode: 'copy',
        pattern: '*.fastq.gz',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/fastq/checksum" },
        pattern: '*.fastq.gz.md5',
        mode: 'copy'
    ]]


    input:
    tuple val(meta), path(reads)

    output:
    path "*.fastq.gz",  includeInputs:true, emit: readsout
    path "*.fastq.gz.md5"                 , emit: md5
    path "versions.yml"                   , emit: versions


    script:
    """

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
    touch tmp.fastq.gz
    touch tmp.fastq.gz.md5
    touch versions.yml
    """
}
