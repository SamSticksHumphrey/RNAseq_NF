
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.


process RSEQC_JUNCTIONSATURATION {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.{pdf,r}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.{pdf.md5,r.md5}'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript
    path "*.md5"                  , emit: md5
    path  "versions.yml"          , emit: versions

    script:
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o ${meta.sample_id}

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_saturation.py --version | sed -e "s/junction_saturation.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.pdf
    touch tmp.r
    touch versions.yml
    """

}
