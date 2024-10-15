process SUBREAD_FEATURECOUNTS {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/featurecounts" },
        mode: 'copy',
        pattern: '*.{txt,summary}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/featurecounts/checksum" },
        mode: 'copy',
        pattern: '*.{txt.md5,summary.md5}'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand), path(annotation)

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "*.md5"                                       , emit: md5
    path "versions.yml"                                , emit: versions

    script:

    def strandedness = 0
    if (strand == 'forward') {
        strandedness = 1
    } else if (strand == 'reverse') {
        strandedness = 2
    }

    """
    featureCounts \\
        -B -C -p \\
        -g gene_biotype \\
        -t exon \\
        -T $task.cpus \\
        -a $annotation \\
        -s $strandedness \\
        -o ${meta.sample_id}_geneBiotype_featureCounts.txt \\
        ${bam}

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$( echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.featureCounts.txt
    touch tmp.featureCounts.txt.summary
    touch versions.yml
    """
}
