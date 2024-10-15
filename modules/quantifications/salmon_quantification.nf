process SALMON_QUANT {
    tag "$meta.sample_id"
    label "process_medium"
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/salmon" },
        mode: 'copy',
        pattern: '*_salmon'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/salmon/checksum" },
        mode: 'copy',
        pattern: '*.{gz.md5,tsv.md5,json.md5,txt.md5,log.md5,sf.md5}'
    ]]

    input:
    tuple val(meta), path(bam), val(strand)
    path  transcript_fasta
    path  gtf


    output:
    tuple val(meta), path("${meta.sample_id}_salmon/") , emit: results
    path "*.md5"                          , emit: md5
    path  "versions.yml"                , emit: versions

    script:
    // May need to use the output of infer experiment to actually infer the libtype if auto is a problem
    """
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        --libType='A' \\
        -t $transcript_fasta \\
        -a $bam \\
        -o ${meta.sample_id}_salmon

    find . -maxdepth 3 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.sample_id}
    touch ${meta.sample_id}/salmon.txt
    touch versions.yml
    """
}
