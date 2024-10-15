#!/bin/bash

# ------------------------------------------------------------------------------
# BuildKallistoIndex.sh  S. Humphrey Jan 2022
# Script to build the indices for the Kallisto pseudoaligner
# ------------------------------------------------------------------------------

# Install Kallisto properly, seems to not work with conda
# 

dataDir=/home/samh/resources/
transcript_fasta=${dataDir}Hsapiens_GRCh38/Homo_sapiens.GRCh38.cdna.all.fa

mkdir -p ${dataDir}kallistoIndex/

kallisto index --make-unique --index ${dataDir}kallistoIndex/kallistoIndex_GRCh38v105.idx \
                $transcript_fasta

