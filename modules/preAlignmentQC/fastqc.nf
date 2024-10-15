process FASTQC {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/fastqc" },
        mode: 'copy',
        pattern: '*.{html,zip}',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/fastqc/checksum" },
        mode: 'copy',
        pattern: '*.{html.md5,zip.md5}'
    ]]


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "*.md5"                   , emit: md5
    path  "versions.yml"           , emit: versions

    script:

    """
    [ ! -f ${meta.sample_id}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.sample_id}_raw_1.fastq.gz
    [ ! -f ${meta.sample_id}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.sample_id}_raw_2.fastq.gz

    fastqc --quiet --threads $task.cpus ${meta.sample_id}_raw_1.fastq.gz ${meta.sample_id}_raw_2.fastq.gz

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    """
    touch fastqc_tmp.html
    touch fastqc_tmp.zip
    touch versions.yml
    """

}
