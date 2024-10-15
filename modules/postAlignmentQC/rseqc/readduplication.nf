
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.

process RSEQC_READDUPLICATION {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.{pdf,xls,r}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.{pdf.md5,xls.md5,r.md5}'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)

    output:
    tuple val(meta), path("*seq.DupRate.xls"), emit: seq_xls
    tuple val(meta), path("*pos.DupRate.xls"), emit: pos_xls
    tuple val(meta), path("*.pdf")           , emit: pdf
    tuple val(meta), path("*.r")             , emit: rscript
    path "*.md5"                          , emit: md5
    path  "versions.yml"                     , emit: versions

    script:
    """
    read_duplication.py \\
        -i $bam \\
        -o ${meta.sample_id}

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_duplication.py --version | sed -e "s/read_duplication.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp_seq.DupRate.xls
    touch tmp_pos.DupRate.xls
    touch tmp.pdf
    touch tmp.r
    touch versions.yml
    """


}
