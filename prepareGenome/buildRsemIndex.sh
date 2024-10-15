#!/bin/bash

# ------------------------------------------------------------------------------
# BuildRsemIndex.sh  S. Humphrey Jan 2022
# Script to build the indices for the Rsem quantification using bowtie aligner
# ------------------------------------------------------------------------------

dataDir=/home/samh/resources/
fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.105.chr.gtf
cpus=30

mkdir -p ${dataDir}rsemIndex/

rsem-prepare-reference --gtf $gtf \
                --num-threads $cpus \
                --bowtie \
                $fasta \
                ${dataDir}rsemIndex/Hsapiens