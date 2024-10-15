process MULTIQC_INDIVIDUAL {
    tag "$meta.sample_id"
    label 'process_low'
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/multiqc" },
        mode: 'copy',
        pattern: '*{_data,_plots,.html}'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/multiqc/checksum" },
        mode: 'copy',
        pattern: '*.md5'
    ]]

    input:
    path multiqc_config
    path software_versions
    tuple val(meta), path(star_log_final), path(fastqc_zip), path(trim_zip), path(trim_log), path(markduplicates_metrics), path(rnaseq_metrics), path(samtools_stats), path(samtools_flagstats), path(samtools_idxstats), path(rsem_stat), path(salmon_results), path(kallisto_results), path(featurecounts_gene_exon), path(dupradar_multiqc), path(featurecounts_multiqc), path(qualimap_multiqc), path(rseqc_bamstat_multiqc), path(rseqc_infer_experiment), path(rseqc_inner_distance), path(rseqc_junction_annotation), path(rseqc_junction_saturation), path(rseqc_read_distribution), path(rseqc_read_duplication), path(rseqc_tin)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*.md5"               , emit: md5
    path "*_plots"             , optional:true, emit: plots

    script:

    title = "\"${params.multiqc_title} : ${meta.sample_id}\""

    """
    multiqc \\
        -f \\
        --title $title \\
        .

    find . -maxdepth 2 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done
    """

    stub:
    """
    touch tmp_multiqc_report.html
    touch tmp_data
    touch tmp_plots
    touch versions.yml
    """
}
