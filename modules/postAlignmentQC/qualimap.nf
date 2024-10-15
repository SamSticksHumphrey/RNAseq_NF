process QUALIMAP_RNASEQ {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/" },
        mode: 'copy',
        pattern: 'qualimap_*',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/qualimap_${meta.sample_id}/checksum" },
        mode: 'copy',
        pattern: '*.{txt.md5,html.md5,png.md5}'
    ]]

    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path  gtf

    output:
    tuple val(meta), path("qualimap_${meta.sample_id}/"), emit: results
    path "*.md5"                       , emit: md5
    path  "versions.yml"               , emit: versions

    script:

    def memory     = task.memory.toGiga() + "G"
    def strandedness = 'non-strand-specific'

// md5 sums for the files with spaces doesn't work. Just running md5sums for the main output files for now

    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    qualimap \\
        --java-mem-size=$memory \\
        rnaseq \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        -pe \\
        -outdir qualimap_${meta.sample_id}

    find . -maxdepth 2 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.sample_id}
    touch ${meta.sample_id}/qualimap.txt
    touch versions.yml
    """
}
