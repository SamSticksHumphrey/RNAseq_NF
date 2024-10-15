process SAMTOOLS_INDEX {
    tag "$meta.sample_id"
    label 'process_low'
    publishDir = [
        path: { "${params.outdir}/${meta.cellline_id}/${meta.sample_id}/star" },
        mode: 'copy',
        pattern: '*.{bai}'
    ]

    input:
    tuple val(meta), path(reads)
    path  fastas

    output:
??

    script:
    sample_id = $meta.sample_id

    """
    bowtie \\
      -1 ${reads[0]}  \\
      -2 ${reads[1]}  \\
      ${viruses_bacteria_index}/viruses_bacteria \\
      --un unmapped.fastq \\
      --phred33-quals \\
      -k 1 \\
      --chunkmbs 128 \\
      --seedmms 0 --seedlen 60 \\
      --seed 123456 \\
      --suppress 1,2,4,5,6,7,8 \\
      --max /dev/null \\
      -p ${task.cpus} \\
      2>> ${sample_id}.viruses_bacteria.align_stats.txt \\
      > aligned_to_genomes.out

    cat aligned_to_genomes.out | \\
      LC_ALL=C \\
      sort \\
      --buffer-size ${task.memory.toGiga()}g \\
      --parallel ${task.cpus} \\
      -T ./tmp/ | \\
      LC_ALL=C \\
      uniq -c | \\
      LC_ALL=C \\
      sort -rn -k 1,1 \\
      -T ./tmp/ |
      sed 's/^ *//g' | sed 's/ /\t/g' \\
      > viruses_bacteria_body.txt

    printf \\
      "Counts_of_mapping_reads\tVirus/Bacteria/Phage\n" \\
      > viruses_bacteria_header.txt

    cat \\
      viruses_bacteria_header.txt \\
      viruses_bacteria_body.txt \\
      > ${sample_id}.viruses_bacteria.aligned_count.txt
    """

    stub:
    """
    touch tmp.bai
    touch versions.yml
    """
}
