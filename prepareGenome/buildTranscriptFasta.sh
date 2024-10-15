#!/bin/bash

# ------------------------------------------------------------------------------
# BuildTranscriptFasta.sh  S. Humphrey Feb 2022
# 
# ------------------------------------------------------------------------------

dataDir=/home/samh/resources/Hsapiens_GRCh38/
genome_fasta=${dataDir}Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=${dataDir}Homo_sapiens.GRCh38.105.chr.gtf



gffread -w ${dataDir}transcripts.fa -g $genome_fasta $gtf
