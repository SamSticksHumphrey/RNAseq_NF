// Run mpileup for NGSCheckMate profileing

process SAMTOOLS_MPILEUP {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
      path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/ngscheckmate/" },
      mode: 'copy',
      pattern: '*.vcf'
      ],[
      path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/qc/ngscheckmate/checksum" },
      mode: 'copy',
      pattern: '*.vcf.md5'
    ]]


    input:
    tuple val(meta), path(bam), path(bai), val(strand)
    path NGSCheckMate_SNPbed
    path fasta

    output:
    path "${outputFile}"     , emit: vcffile
    path "*.md5"             , emit: md5
    path  "versions.yml"     , emit: versions

    script:
    outputFile = "${meta.sample_id}_NGSCheckMate.vcf"

    """
    bcftools mpileup \\
        -f $fasta  \\
        -R $NGSCheckMate_SNPbed  \\
        --max-depth 500  \\
        -Ou  \\
        -I ${bam} \\
          | bcftools call -c -Ov -o $outputFile

    find . -maxdepth 1 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | awk -F 'Copyright' '{print \$1}')
    END_VERSIONS
    """

    stub:
    outputFile = "${meta.sample_id}_RNA.vcf"
    """
    touch $outputFile
    touch versions.yml
    """
}
