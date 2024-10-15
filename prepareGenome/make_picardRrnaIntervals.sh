#!/usr/bin/env bash

# make_rRNA.sh
# Kamil Slowikowski
# December 12, 2014
# Modified: Arindam Ghosh (July 24, 2019 )
# Modified: Sam Humphrey (March, 2022 )
#
# Referenc Genome: GRCh38 Ensembl release 105
#
# 1. Prepare chromosome sizes file from fasta sequence if needed.
# 2. Make an interval_list file suitable for CollectRnaSeqMetrics.jar.

resources_dir=/home/samh/resources/Hsapiens_GRCh38/

chrom_sizes=${resources_dir}Homo_sapiens.GRCh38.dna.primary_assembly.fa.sizes

# rRNA interval_list file -------------------------------------------------

# Genes from Ensembl.

genes=${resources_dir}Homo_sapiens.GRCh38.105.chr.gtf

# Output file suitable for Picard CollectRnaSeqMetrics.jar.

rRNA=${resources_dir}Homo_sapiens.GRCh38.105.rRNA.interval_list

# Sequence names and lengths. (Must be tab-delimited.)
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' $chrom_sizes | \
    grep -v _ \
>> $rRNA

# Intervals for rRNA transcripts.
egrep -w 'gene_biotype "rRNA"|gene_biotype "rRNA_pseudogene"|gene_biotype "Mt_rRNA"' $genes | \
    awk '$3 == "gene"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /gene_id "([^"]+)"/ or die "no gene_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA
