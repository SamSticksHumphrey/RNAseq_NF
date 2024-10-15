#!/usr/bin/env Rscript


# -----------------------------------------------------------------------------------------------------------------------
#
# We could alternatively generate counts from abundances, using the argument countsFromAbundance, scaled to library size, "scaledTPM",
#    or additionally scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM".
# Using either of these approaches, the counts are not correlated with length, and so the length matrix should not be provided as an offset for downstream analysis packages.
# As of tximport version 1.10, we have added a new countsFromAbundance option "dtuScaledTPM".
# This scaling option is designed for use with txOut=TRUE for differential transcript usage analyses.
# See ?tximport for details on the various countsFromAbundance options.
#
# -----------------------------------------------------------------------------------------------------------------------

#library(rtracklayer)
#library(SummarizedExperiment)
library(tximport)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Usage: salmon_tximport.r <coldata> <salmon_out>", call.=FALSE)
}

path = args[1]
gtf_file  = args[2]
prefix = args[3]


# Setup files
gtf <- rtracklayer::import(gtf_file)
tx_rowdata <- gtf %>% as_tibble() %>% filter(type == 'transcript') %>% select( gene_id, gene_name, gene_biotype, transcript_id, transcript_name, transcript_biotype ) %>% distinct()
tx2gene <- tx_rowdata %>% select(transcript_id, gene_id)
gene_rowdata <- tx_rowdata %>% select(-starts_with('transcript')) %>% distinct()

salmon_output <- list.files(path, pattern = "quant.sf", recursive = T, full.names = T)

# run tx import
txi_gene <- tximport::tximport(salmon_output, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
txi_gene_scaledTPM <- tximport::tximport(salmon_output, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "scaledTPM")
txi_gene_lengthScaledTPM <- tximport::tximport(salmon_output, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

txi_tx = tximport::tximport(salmon_output, type = "salmon", txOut = TRUE)
txi_dtuScaledTPM <- tximport::tximport(salmon_output, type = "salmon", txOut = TRUE, tx2gene = tx2gene, countsFromAbundance = "dtuScaledTPM")

# Abundance is the TPM table from the original salmon quantification, but gene-level summarized
gene_tpm <- as_tibble(txi_gene$abundance, rownames = 'gene_id') %>%
    full_join(gene_rowdata, by = 'gene_id') %>%
    select(starts_with('gene'), everything()) %>% 
    dplyr::rename(!!prefix := V1) 

write.table(gene_tpm , paste0(prefix, "_gene_tpms.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

# The count table from the original salmon quantification, but gene-level summarized.
gene_counts <- as_tibble(txi_gene$counts, rownames = 'gene_id') %>%
    full_join(gene_rowdata, by = 'gene_id') %>%
    select(starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1)

write.table(gene_counts, paste0(prefix, "_gene_counts.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

# Gene lengths from the original salmon quantification
gene_lengths <- as_tibble(txi_gene$length, rownames = 'gene_id') %>%
    full_join(gene_rowdata, by = 'gene_id') %>%
    select(starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(gene_lengths, paste0(prefix, "_gene_lengths.tsv"), sep="\t", quote=FALSE, row.names = FALSE)


# Counts will be generated by using the TPM value * featureLength * library size.
gene_lengthScaledTPM <- as_tibble(txi_gene_lengthScaledTPM$counts, rownames = 'gene_id') %>%
    full_join(gene_rowdata, by = 'gene_id') %>%
    select(starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(gene_lengthScaledTPM, paste0(prefix, "_gene_counts_scaledByLibrarySizeAndGeneLength.tsv"),
            sep="\t", quote=FALSE, row.names = FALSE)



# The count table from the original salmon quantification at transcript level.
tx_counts <- as_tibble(txi_tx$counts, rownames = 'transcript_id') %>%
    # full_join(tx_rowdata, by = 'transcript_id') %>%
    select(starts_with('transcript'), starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(tx_counts, paste0(prefix, "_transcript_counts.tsv"), sep="\t", quote=FALSE, row.names = FALSE)


# The abundance from the original salmon quantification at transcript level.
tx_TPM <- as_tibble(txi_tx$abundance, rownames = 'transcript_id') %>%
    #  full_join(tx_rowdata, by = 'transcript_id') %>%
    select(starts_with('transcript'), starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(tx_TPM, paste0(prefix, "_transcript_tpms.tsv"), sep="\t", quote=FALSE, row.names = FALSE)

# The transcript lengths from the original salmon quantification.
tx_lengths <- as_tibble(txi_tx$length, rownames = 'transcript_id') %>%
    # full_join(tx_rowdata, by = 'transcript_id') %>%
    select(starts_with('transcript'), starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(tx_lengths, paste0(prefix, "_transcript_lengths.tsv"), sep="\t", quote=FALSE, row.names = FALSE)


# This scaling option is designed for use with txOut=TRUE for differential transcript usage analyses.
txi_dtuScaledTPM <- as_tibble(txi_dtuScaledTPM$counts, rownames = 'transcript_id') %>%
    #   full_join(tx_rowdata, by = 'transcript_id') %>%
    select(starts_with('transcript'), starts_with('gene'), everything())%>% 
    dplyr::rename(!!prefix := V1) 

write.table(tx_lengths, paste0(prefix, "_transcript_scaledForDTU.tsv"), sep="\t", quote=FALSE, row.names = FALSE)


# Print sessioninfo to standard out
#citation(rtracklayer)
#citation(tximport)

sessionInfo()
