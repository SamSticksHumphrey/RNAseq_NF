def VERSION = '377' // Version information not provided by tool on CLI

//Convert a bedGraph file to bigWig format.
process BEDGRAPH_TO_BIGWIG {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig" },
            mode: 'copy',
            pattern: '*.bigWig'
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig/checksum" },
            mode: 'copy',
            pattern: '*.bigWig.md5'
        ]]

    input:
    tuple val(meta), path(bedgraph_forward), path(bedgraph_reverse)
    path  sizes

    output:
    tuple val(meta), path("*forward.bigWig"), path("*reverse.bigWig"), emit: bigwigs
    path "*.md5"             , emit: md5
    path "versions.yml"                                             , emit: versions

    script:
    """
    bedGraphToBigWig \\
        $bedgraph_forward \\
        $sizes \\
        ${meta.sample_id}_forward.bigWig

    bedGraphToBigWig \\
        $bedgraph_reverse \\
        $sizes \\
        ${meta.sample_id}_reverse.bigWig


    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """

    stub:
    """
    touch tmp_forward.bigWig
    touch tmp_reverse.bigWig
    touch versions.yml
    """
}
