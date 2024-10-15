process RSEM_CALCULATEEXPRESSION {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/rsem" },
            mode: 'copy',
            pattern: "*.{stat,results}"
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/rsem" },
            mode: 'copy',
            pattern: "*.bam",
            enabled: params.save_all_intermediates
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/rsem/checksum" },
            mode: 'copy',
            pattern: "*.{stat.md5,results.md5,theta.md5,model.md5,cnt.md5}"
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/rsem/checksum" },
            mode: 'copy',
            pattern: "*.bam.md5",
            enabled: params.save_all_intermediates
        ]]



    input:
    tuple val(meta), path(bam), val(strand)
    path index

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    path "*.md5"                               , emit: md5
    path  "versions.yml"                       , emit: versions

    script:
    sample_id   =  "${meta.sample_id}"

    """
    rsem-calculate-expression --alignments \\
      -p $task.cpus \\
      --paired-end \\
      --strandedness $strand \\
      --no-bam-output \\
      ${bam} \\
      ${index}/Hsapiens \\
      ${sample_id}

    find . -maxdepth 2 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    sample_id   =  "${meta.sample_id}"
    """
    touch ${sample_id}.genes.results
    touch ${sample_id}.isoforms.results
    touch ${sample_id}.stat
    touch ${sample_id}.log
    touch versions.yml
    """
}
