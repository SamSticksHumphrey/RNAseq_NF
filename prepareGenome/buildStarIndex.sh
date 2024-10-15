#!/bin/bash

# ------------------------------------------------------------------------------
# BuildStarIndex.sh  S. Humphrey Jan 2022
# Script to build the indices for the Star aligner
# ------------------------------------------------------------------------------

dataDir=/home/samh/resources/
fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.105.chr.gtf
cpus=30

mkdir ${dataDir}starIndex

STAR \
    --runMode genomeGenerate \
    --genomeDir ${dataDir}starIndex/ \
    --genomeFastaFiles $fasta \
    --sjdbGTFfile $gtf \
    --runThreadN $cpus \
    --sjdbOverhang 160
