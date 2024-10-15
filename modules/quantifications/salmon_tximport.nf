process SALMON_TXIMPORT {
    tag "$meta.sample_id"
    label "process_medium"
    publishDir = [[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/salmon" },
        mode: 'copy',
        pattern: '*.tsv'
    ],[
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/salmon/checksum" },
        mode: 'copy',
        pattern: '*.tsv.md5'
    ]]

    input:
    tuple val(meta), path("${meta.sample_id}/*")
    path gtf

    output:
    path "${meta.sample_id}_gene_tpms.tsv"                         , emit: gene_tpms
    path "${meta.sample_id}_gene_counts.tsv"                       , emit: gene_counts
    path "${meta.sample_id}_gene_lengths.tsv"                      , emit: gene_lengths
    path "${meta.sample_id}_gene_counts_scaledByLibrarySizeAndGeneLength.tsv"   , emit: gene_counts_scaledByLibrarySizeAndGeneLength

    path "${meta.sample_id}_transcript_counts.tsv"     , emit: transcript_counts
    path "${meta.sample_id}_transcript_tpms.tsv"       , emit: transcript_tpms
    path "${meta.sample_id}_transcript_lengths.tsv"    , emit: transcript_lengths
    path "${meta.sample_id}_transcript_scaledForDTU.tsv"   , emit: transcript_scaledForDTU
    path "*.md5"                          , emit: md5
    path "versions.yml"                          , emit: versions

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    salmon_tximport.r ${meta.sample_id} $gtf ${meta.sample_id}

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.sample_id}_gene_tpms.tsv"
    touch "${meta.sample_id}_gene_counts.tsv"
    touch "${meta.sample_id}_gene_lengths.tsv"
    touch "${meta.sample_id}_gene_counts_scaledByLibrarySizeAndGeneLength.tsv"

    touch "${meta.sample_id}_transcript_counts.tsv"
    touch "${meta.sample_id}_transcript_tpms.tsv"
    touch "${meta.sample_id}_transcript_lengths.tsv"
    touch "${meta.sample_id}_transcript_scaledForDTU.tsv"
    touch "versions.yml"
    """
}
