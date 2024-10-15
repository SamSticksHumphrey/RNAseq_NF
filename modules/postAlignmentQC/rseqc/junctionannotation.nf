
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.

process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.{pdf,bed,xls,log}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.{pdf.md5,Interact.bed.md5,junction.bed.md5,xls.md5,log.md5}'
    ]]


    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  bed

    output:
    tuple val(meta), path("*.xls")         , emit: xls
    tuple val(meta), path("*.r")           , emit: rscript
    tuple val(meta), path("*.log")         , emit: log
    tuple val(meta), path("*.junction.bed"), emit: bed
    tuple val(meta), path("*.Interact.bed"), emit: interact_bed
    tuple val(meta), path("*junction.pdf") , emit: pdf
    tuple val(meta), path("*events.pdf")   , emit: events_pdf
    path "*.md5"                           , emit: md5
    path  "versions.yml"                   , emit: versions

    script:
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o ${meta.sample_id} \\
        2> ${meta.sample_id}.junction_annotation.log

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_annotation.py --version | sed -e "s/junction_annotation.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch jctAnn.xls
    touch jctAnn.r
    touch jctAnn.log
    touch jctAnn.junction.bed
    touch jctAnn.Interact.bed
    touch jctAnn_junction.pdf
    touch jctAnn_events.pdf
    touch versions.yml
    """

}
