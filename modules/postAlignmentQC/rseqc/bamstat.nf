
// RSeQC package provides a number of useful modules that can comprehensively evaluate
// high throughput sequence data especially RNA-seq data.

process RSEQC_BAMSTAT {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.bam_stat.txt'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.bam_stat.txt.md5'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt
    path "*.md5"                          , emit: md5
    path  "versions.yml"                   , emit: versions

    script:
    """
    bam_stat.py \\
        -i $bam \\
        > ${meta.sample_id}.bam_stat.txt
    
    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(bam_stat.py --version | sed -e "s/bam_stat.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.bam_stat.txt
    touch versions.yml
    """

}
