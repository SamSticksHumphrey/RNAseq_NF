process UNTAR {
    tag "$archive"
 //   label 'process_low'

    input:
    path archive

    output:
    path "$untar"      , emit: untar
    path "versions.yml", emit: versions

    script:
    untar        = archive.toString() - '.tar.gz'
    """
    tar \\
        -xzvf \\
        $archive \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        untar: \$(echo \$(tar --version 2>&1) | sed 's/^.*(GNU tar) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch $untar
    touch versions.yml
    """

}
