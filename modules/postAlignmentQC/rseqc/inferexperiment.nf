process RSEQC_INFEREXPERIMENT {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.infer_experiment.txt'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.infer_experiment.txt.md5'
    ]]

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: inferexperiment_txt
    path "*.md5"                                   , emit: md5
    path  "versions.yml"                           , emit: versions

    script:

    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        > ${meta.sample_id}.infer_experiment.txt

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.txt
    touch versions.yml
    """
}
