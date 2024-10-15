process SAMTOOLS_SORT_INDEX {
    tag "$meta.sample_id"
    label 'process_high'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam_bai
    path "*.md5"                   , emit: md5
    path  "versions.yml"          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    samtools sort -@ $task.cpus -o ${prefix}.bam -T ${prefix} $bam

    samtools index -@ ${task.cpus-1} ${prefix}.bam

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
    touch tmp_sorted.bam
    touch versions.yml
    """
}
