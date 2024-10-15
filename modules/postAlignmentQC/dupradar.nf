process DUPRADAR {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/dupradar" },
            mode: 'copy',
            pattern: "*.{pdf,txt}"
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/dupradar/checksum" },
            mode: 'copy',
            pattern: "*.{pdf.md5,txt.md5}"
        ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  gtf

    output:
    tuple val(meta), path("*.pdf")    , emit: pdf
    tuple val(meta), path("*.txt")    , emit: txt
    tuple val(meta), path("*_mqc.txt"), emit: multiqc
    path "*.md5"                      , emit: md5
    path "versions.yml"               , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/

    def strandedness = 0
    if (strand == 'forward') {
        strandedness = 1
    } else if (strand == 'reverse') {
        strandedness = 2
    }

    """
    dupradar.r \\
        $bam \\
        ${meta.sample_id} \\
        $gtf \\
        $strandedness \\
        'paired' \\
        $task.cpus

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-dupradar: \$(Rscript -e "library(dupRadar); cat(as.character(packageVersion('dupRadar')))")
    END_VERSIONS
    """

    stub:
    """
    touch tmp.pdf
    touch tmp.txt
    touch tmp_mqc.txt
    touch versions.yml
    """
}
