process MULTIQC_BATCH {
    label 'process_low'
    publishDir = [
        path: { "${params.batch}/multiqc" },
        mode: 'copy'
        ]

    input:
    path multiqc_config
    path software_versions
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , emit: plots

    script:

    title = "\"${params.multiqc_title} : ${params.jobname}\""

    """
    multiqc \\
        -f \\
        --title $title \\
        .
    """

    stub:
    """
    touch tmp_multiqc_report.html
    touch tmp_data
    touch tmp_plots
    touch versions.yml
    """
}
