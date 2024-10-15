process MULTIQC_CUSTOM_BIOTYPE {
    tag "$meta.sample_id"
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/featurecounts" },
        mode: 'copy',
        pattern: '*.tsv'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/featurecounts/checksum" },
        mode: 'copy',
        pattern: '*.tsv.md5'
    ]]

    input:
    tuple val(meta), path(count)
    path  header

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "*.md5"             , emit: md5
    path "versions.yml"           , emit: versions

    script:

    """
    cut -f 1,7 $count | tail -n +3 | cat $header - >> ${meta.sample_id}.biotype_counts_mqc.tsv

    mqc_features_stat.py \\
        ${meta.sample_id}.biotype_counts_mqc.tsv \\
        -s $meta.sample_id \\
        -f rRNA \\
        -o ${meta.sample_id}.biotype_counts_rrna_mqc.tsv

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch tmp.tsv
    touch versions.yml
    """

}
