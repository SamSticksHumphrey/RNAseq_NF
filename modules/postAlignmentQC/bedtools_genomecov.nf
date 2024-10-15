process BEDTOOLS_GENOMECOV {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig/genomecov" },
        enabled: params.save_all_intermediates,
        mode: 'copy',
        pattern: '*.bedGraph'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/bedGraph_bigwig/genomecov" },
        enabled: params.save_all_intermediates,
        mode: 'copy',
        pattern: '*.bedGraph.md5'
    ]]


    input:
    tuple val(meta), path(bam), path(bai), val(strand)

    output:
    tuple val(meta), path("*_forward.bedGraph"), path("*_reverse.bedGraph"), emit: bedgraphs
    path "*.md5"                                        , emit: md5
    path "versions.yml"                        , emit: versions

    script:

    def sample_id = "${meta.sample_id}"

    def sample_id_forward = "${sample_id}_forward"
    def sample_id_reverse = "${sample_id}_reverse"
    if (strand == 'reverse') {
        sample_id_forward = "${sample_id}_reverse"
        sample_id_reverse = "${sample_id}_forward"
      }

    """
    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand + \\
        -split -du \\
        | bedtools sort > ${sample_id_forward}.bedGraph

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -strand - \\
        -split -du \\
        | bedtools sort > ${sample_id_reverse}.bedGraph

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp_forward.bedGraph
    touch tmp_reverse.bedGraph
    touch versions.yml
    """

}
