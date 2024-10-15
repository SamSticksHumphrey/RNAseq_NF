process TRIMGALORE {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/trimgalore" },
            mode: 'copy',
            pattern: "*.{html,zip,txt}"
        ],  [
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/trimgalore" },
            mode: 'copy',
            pattern: "*.fq.gz",
            enabled: params.save_all_intermediates
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/trimgalore/checksum" },
            mode: 'copy',
            pattern: "*.{html.md5,zip.md5,txt.md5}"
        ], [
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/trimgalore/checksum" },
            mode: 'copy',
            pattern: "*.fq.gz.md5",
            enabled: params.save_all_intermediates
        ]]

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")    , emit: reads
    tuple val(meta), path("*report.txt"), emit: log
    path "versions.yml"                 , emit: versions

    tuple val(meta), path("*.html"), emit: html
    path "*.md5"                          , emit: md5
    tuple val(meta), path("*.zip") , emit: zip

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Clipping presets have to be evaluated in the context of SE/PE
  //  def c_r1   = params.clip_r1 > 0             ? "--clip_r1 ${params.clip_r1}"                         : ''
  //  def c_r2   = params.clip_r2 > 0             ? "--clip_r2 ${params.clip_r2}"                         : ''
  //  def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
  //  def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
  // Add to trim_galore if needed:
  // $c_r1 \\
  // $c_r2 \\
  // $tpc_r1 \\
  // $tpc_r2 \\


    // Added soft-links to original fastqs for consistent naming in MultiQC

    """
    [ ! -f  ${meta.sample_id}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.sample_id}_1.fastq.gz
    [ ! -f  ${meta.sample_id}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.sample_id}_2.fastq.gz

    trim_galore \\
        --cores $cores \\
        --paired \\
        --gzip \\
        --fastqc \\
        --clip_R1 10 \\
        --clip_R2 10 \\
        ${meta.sample_id}_1.fastq.gz \\
        ${meta.sample_id}_2.fastq.gz


    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS


    """

    stub:
    """
    touch trim_tmp.fq.gz
    touch trim_tmp.report.txt
    touch trim_tmp.html
    touch trim_tmp.zip
    touch versions.yml
    """
}
