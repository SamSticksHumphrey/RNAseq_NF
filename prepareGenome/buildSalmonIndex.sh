#!/bin/bash

# ------------------------------------------------------------------------------
# BuildSalmnonIndex.sh  S. Humphrey Jan 2022
# Script to build the gentrome and indices for Salmon aligner
# ------------------------------------------------------------------------------

dataDir=/home/samh/resources/
genome_fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
transcript_fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.cdna.all.fa
nctranscript_fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.ncrna.fa
gtf=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.105.chr.gtf
gentrome_fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.gentrome.fa
cpus=30

mkdir -p ${dataDir}salmonIndex

grep '^>' $genome_fasta | cut -d ' ' -f 1 > ${dataDir}Hsapiens_GRCh38/decoys.txt
sed -i.bak -e 's/>//g' ${dataDir}Hsapiens_GRCh38/decoys.txt

cat $transcript_fasta $nctranscript_fasta $genome_fasta > $gentrome_fasta

salmon index \
    --threads $cpus \
    -t $gentrome_fasta \
    -d ${dataDir}Hsapiens_GRCh38/decoys.txt \
    -i ${dataDir}salmonIndex
