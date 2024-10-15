process STAR_ALIGN {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/star" },
            mode: 'copy',
            pattern: '*.{out,tab,bam}'
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/star/checksum" },
            mode: 'copy',
            pattern: '*.{out.md5,tab.md5,bam.md5}'
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/star/unmapped" },
            mode: 'copy',
            pattern: '*.fastq.gz',
            enabled: params.save_all_intermediates
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/star/unmapped/checksum" },
            mode: 'copy',
            pattern: '*.fastq.gz.md5',
            enabled: params.save_all_intermediates
        ]]

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*Aligned.out.bam')           , emit: bam
    tuple val(meta), path('*Log.final.out')             , emit: log_final
    tuple val(meta), path('*Log.out')                   , emit: log_out
    tuple val(meta), path('*Log.progress.out')          , emit: log_progress
    tuple val(meta), path('*toTranscriptome.out.bam')   , emit: bam_transcript
    tuple val(meta), path('*.tab')                      , emit: tab
    path "versions.yml"                                 , emit: versions
    path "*.md5"                                        , emit: md5

    script:
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${meta.sample_id}_ \\
        --quantMode TranscriptomeSAM \\
        --outSAMtype BAM Unsorted \\
        --alignSJDBoverhangMin 1 \\
        --outSAMmultNmax 10 \\
        --outFilterMultimapNmax 10 \\
        --twopassMode Basic \\
        --outSAMattributes NH HI AS NM MD


    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    """
    touch tmp_Aligned.out.bam
    touch tmp_Log.final.out
    touch tmp_Log.out
    touch tmp_Log.progress.out
    touch versions.yml

    touch tmp_toTranscriptome.out.bam
    touch tmp_Aligned.unsort.out.bam
    touch tmp_fastq.gz
    touch tmp_.tab
    """
}
