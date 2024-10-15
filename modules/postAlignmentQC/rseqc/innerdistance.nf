
//  RSeQC package provides a number of useful modules that can comprehensively evaluate
//    high throughput sequence data especially RNA-seq data.

process RSEQC_INNERDISTANCE {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc" },
        mode: 'copy',
        pattern: '*.{txt,pdf,r}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/rseqc/checksum" },
        mode: 'copy',
        pattern: '*.{txt.md5,pdf.md5,r.md5}'
    ]]


    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  bed

    output:
    tuple val(meta), path("*distance.txt"), emit: distance
    tuple val(meta), path("*freq.txt")    , emit: freq
    tuple val(meta), path("*mean.txt")    , emit: mean
    tuple val(meta), path("*.pdf")        , emit: pdf
    tuple val(meta), path("*.r")          , emit: rscript
    path "*.md5"                          , emit: md5
    path  "versions.yml"                  , emit: versions

    script:
    """
    inner_distance.py \\
        -i $bam \\
        -r $bed \\
        -o ${meta.sample_id} \\
        > stdout.txt
    head -n 2 stdout.txt > ${meta.sample_id}.inner_distance_mean.txt

    rm stdout.txt

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
    END_VERSIONS
    """


    stub:
    """
    touch tmp_distance.txt
    touch tmp_freq.txt
    touch tmp_mean.txt
    touch tmp.pdf
    touch tmp.r
    touch versions.yml
    """
}
