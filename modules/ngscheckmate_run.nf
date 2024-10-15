// Overall NGSCheckMate profiling
// this will need a filter over the intervals for targeted sequencing

process RUN_NGSCHECKMATE {
    label 'process_low'
    container '449547545634.dkr.ecr.eu-west-2.amazonaws.com/ngscheckmate:1.0'
    publishDir = [
        path: { "${params.batch}/ngscheckmate/" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename.endsWith(".md5") ? "checksum/$filename" : filename }
    ]

    input:
    path multiqc_files
    path snp_bed

    output:
    path "*"                , emit: output
    path "*.md5"            , emit: md5
    path  "versions.yml"    , emit: versions

    script:
    def NGSCHECKMATE_VERSION = '1.0'

    """

    # Run NGSCheckMate - using a local version of the tool
    python2.7 \${NGSCHECKMATE}ncm.py \
                    -V \
                    -d . \
                    -bed $snp_bed \
                    -O .


    find . -maxdepth 2 -type f | \
    while read FILE
    do
        TMPB=\$(basename \$FILE)
        md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NGSCheckMate: $NGSCHECKMATE_VERSION
    END_VERSIONS
    """

    stub:
    """
        touch versions.yml
        touch tmp.md5
    """
}
