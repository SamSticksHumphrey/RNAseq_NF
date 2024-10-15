//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA } from '../modules/prepareGenome/gunzip'
include { GUNZIP as GUNZIP_GTF } from '../modules/prepareGenome/gunzip'
include { GUNZIP as GUNZIP_TX } from '../modules/prepareGenome/gunzip'
include { UNTAR } from '../modules/prepareGenome/untar'
include { GTF2BED              } from '../modules/prepareGenome/gtf2bed'
include { GTF_GENE_FILTER      } from '../modules/prepareGenome/gtf_gene_filter'
include { GET_CHROM_SIZES      } from '../modules/prepareGenome/get_chrom_sizes'

/*
include { SALMON_INDEX } from '../modules/prepareGenome/buildSalmonIndex'
include { RSEM_INDEX } from '../modules/prepareGenome/buildRsemIndex'
include { STAR_INDEX } from '../modules/prepareGenome/buildStarIndex'
include { KALLISTO_INDEX } from '../modules/prepareGenome/buildKallistoIndex'
*/

workflow PREPARE_GENOME {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.genome_fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA( [ [:], params.genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.genome_fasta)
    }


    //
    // Uncompress GTF annotation file
    //
    if (params.gtf.endsWith('.gz')) {
        ch_gtf      = GUNZIP_GTF( [ [:], params.gtf ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = file(params.gtf)
    }

    //
    // create gene bed from GTF
    //
    ch_gene_bed = GTF2BED ( ch_gtf ).bed
    ch_versions = ch_versions.mix(GTF2BED.out.versions)


    //
    // Uncompress transcript fasta file and join with nc-RNA
    //
    if (params.transcript_fasta.endsWith('.gz')) {
        ch_transcript_fasta = GUNZIP_TX( [:], params.transcript_fasta ).gunzip.map { it[1] }
        ch_versions         = ch_versions.mix(GUNZIP_TX.out.versions)
    } else {
        ch_transcript_fasta = file(params.transcript_fasta)
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = GET_CHROM_SIZES ( ch_fasta ).sizes
    ch_versions    = ch_versions.mix(GET_CHROM_SIZES.out.versions)




    // Uncompress STAR index or generate from scratch if required
    ch_star_index = Channel.empty()
    if (params.star_index) {
        if (params.star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR( params.star_index ).untar
            ch_versions   = ch_versions.mix(UNTAR.out.versions)
        } else {
            ch_star_index = file(params.star_index)
        }
    } else {
      println 'Need to include the star index building script'
      //  ch_star_index = STAR_INDEX ( ch_fasta, ch_gtf ).index
      //  ch_versions   = ch_versions.mix(STAR_INDEX.out.versions)
    }

    // Uncompress RSEM index or generate from scratch if required
    ch_rsem_index = Channel.empty()
        if (params.rsem_index) {
            if (params.rsem_index.endsWith('.tar.gz')) {
                ch_rsem_index = UNTAR( params.rsem_index ).untar
                ch_versions   = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_rsem_index = file(params.rsem_index)
            }
        } else {
            println 'Need to include the rsem index building script'
          //  ch_rsem_index = RSEM_INDEX ( ch_fasta, ch_gtf ).index
          //  ch_versions   = ch_versions.mix(RSEM_INDEX.out.versions)
        }


    // Uncompress Kallisto index or generate from scratch if required
    ch_kallisto_index = Channel.empty()
        if (params.kallisto_index) {
            if (params.kallisto_index.endsWith('.tar.gz')) {
                ch_kallisto_index = UNTAR( params.kallisto_index ).untar
                ch_versions   = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_kallisto_index = file(params.kallisto_index)
            }
        } else {
            println 'Need to include the kallisto index building script'
          //  ch_kallisto_index = KALLISTO_INDEX ( ch_fasta, ch_gtf ).index
          //  ch_versions   = ch_versions.mix(KALLISTO_INDEX.out.versions)
        }


    // Uncompress Salmon index or generate from scratch if required
    ch_salmon_index = Channel.empty()
        if (params.salmon_index) {
            if (params.salmon_index.endsWith('.tar.gz')) {
                ch_salmon_index = UNTAR ( params.salmon_index ).untar
                ch_versions     = ch_versions.mix(UNTAR.out.versions)
            } else {
                ch_salmon_index = file(params.salmon_index)
            }
        } else {
            println 'Need to include the kallisto index building script'
        //    ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
        //    ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
    }


    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    gene_bed         = ch_gene_bed         //    path: gene.bed
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    star_index       = ch_star_index       //    path: star/index/
    rsem_index       = ch_rsem_index       //    path: rsem/index/
    kallisto_index   = ch_kallisto_index     //    path: hisat2/index/
    salmon_index     = ch_salmon_index     //    path: salmon/index/

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
