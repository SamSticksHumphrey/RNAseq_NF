process SAMTOOLS_ALLSTATS {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/samtools" },
        mode: 'copy',
        pattern: '*stats'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/samtools/checksum" },
        mode: 'copy',
        pattern: '*stats.md5'
    ]]

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.flagstats"), emit: flagstat
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path  "versions.yml"            , emit: versions
    path "*.md5"                    , emit: md5

    script:

    """
    samtools stats --threads ${task.cpus-1} --ref-seq ${fasta} ${bam} > ${bam}.stats

    samtools idxstats $bam > ${bam}.idxstats

    samtools flagstat --threads ${task.cpus-1} $bam > ${bam}.flagstats

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch tmp.stats
    touch versions.yml
    """
}
