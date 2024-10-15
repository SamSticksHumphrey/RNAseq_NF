process KALLISTO_ALIGN {
    tag "$meta.sample_id"
    label 'process_medium'
    publishDir = [[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/kallisto" },
            mode: 'copy',
            pattern: '*.{h5,tsv,json,log}'
        ],[
            path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/kallisto/checksum" },
            mode: 'copy',
            pattern: '*.{h5.md5,tsv.md5,json.md5,log.md5}'
        ]]

    input:
    tuple val(meta), path(reads), val(strand)
    path  index

    output:
    tuple val(meta), path('abundance.h5')    , emit: abundace_h5
    tuple val(meta), path('abundance.tsv')   , emit: abundace_tsv
    tuple val(meta), path('run_info.json')   , emit: run_info
    tuple val(meta), path('*_kallisto.log')  , emit: log
    path "*.md5"                          , emit: md5
    path "versions.yml"                      , emit: versions


    script:
    """
    [ ! -f ${meta.sample_id}_1.fq.gz ] && ln -s ${reads[0]} ${meta.sample_id}_1.fastq.gz
    [ ! -f ${meta.sample_id}_2.fq.gz ] && ln -s ${reads[1]} ${meta.sample_id}_2.fastq.gz

    kallisto quant \\
        -i $index \\
        -o ./ \\
        -t ${task.cpus} \\
        ${meta.sample_id}_1.fastq.gz ${meta.sample_id}_2.fastq.gz &> ${meta.sample_id}_kallisto.log

    find . -maxdepth 2 -type f | \
    while read FILE
    do
      TMPB=\$(basename \$FILE)
      md5sum \$FILE > \$TMPB.md5
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(kallisto version)
    END_VERSIONS
    """

    stub:
    """
    touch abundance.h5
    touch abundance.tsv
    touch run_info.json
    touch ${meta.sample_id}_kallisto.log
    touch versions.yml
    """

    }
