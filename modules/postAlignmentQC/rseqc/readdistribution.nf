
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.

process RSEQC_READDISTRIBUTION {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.read_distribution.txt'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.read_distribution.txt.md5'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  bed

    output:
    tuple val(meta), path("*.read_distribution.txt"), emit: txt
    path "*.md5"                          , emit: md5
    path  "versions.yml"                            , emit: versions

    script:
    """
    read_distribution.py \\
        -i $bam \\
        -r $bed \\
        > ${meta.sample_id}.read_distribution.txt

    
    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_distribution.py --version | sed -e "s/read_distribution.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.read_distribution.txt
    touch versions.yml
    """
}
