def VERSION = '377' // Version information not provided by tool on CLI

//description: Remove lines from bed file that refer to off-chromosome locations.

process BEDCLIP {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig" },
        pattern: '*.bedGraph',
        mode: 'copy',
        enabled: params.save_all_intermediates
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig/checksum" },
        pattern: '*.bedGraph.md5',
        mode: 'copy',
        enabled: params.save_all_intermediates
    ]]

    input:
    tuple val(meta), path(bedgraph_forward), path(bedgraph_reverse)
    path  sizes

    output:
    tuple val(meta), path("*bedClip_forward.bedGraph"), path("*bedClip_reverse.bedGraph"), emit: bedClip_graphs
    path "*.md5"                                       , emit: md5
    path "versions.yml"                , emit: versions

    script:
    """
    bedClip \\
        $bedgraph_forward \\
        $sizes \\
        ${meta.sample_id}_bedClip_forward.bedGraph

    bedClip \\
        $bedgraph_reverse \\
        $sizes \\
        ${meta.sample_id}_bedClip_reverse.bedGraph

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
    touch tmp2_bedClip_forward.bedGraph
    touch tmp2_bedClip_reverse.bedGraph
    touch versions.yml
    """
}
