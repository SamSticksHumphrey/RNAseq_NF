
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.

process RSEQC_TIN {
    tag "$meta.sample_id"
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.{txt,xls}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.{txt.md5,xls.md5}'
    ]]


    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.txt"), emit: tin_multiqc
    tuple val(meta), path("*.xls"), emit: tin_xls
    path "*.md5"                  , emit: md5
    path "versions.yml"           , emit: versions

    script:
    """
    [ ! -f ${meta.sample_id}.bam ]     && ln -s ${bam} ${meta.sample_id}.bam
    [ ! -f ${meta.sample_id}.bam.bai ] && ln -s ${bai} ${meta.sample_id}.bam.bai

    tin.py \\
        -i ${meta.sample_id}.bam \\
        -r $bed

    sed -i 's/${meta.sample_id}.bam/${meta.sample_id}/g' ${meta.sample_id}.summary.txt

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(tin.py --version | sed -e "s/tin.py //g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.txt
    touch tmp.xls
    touch versions.yml
    """
}
