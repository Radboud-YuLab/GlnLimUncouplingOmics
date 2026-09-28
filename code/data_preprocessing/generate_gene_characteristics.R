## Script that mines gene features
# written by: Cyriel Huijer
# last updated: 19 June 2026


library(tidyverse)
library(biomaRt)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(cubar)
library(Biostrings)
library(ORFik)
library(data.table)
library(Peptides)

# clear environment
rm(list=ls())

# Set working directory to /code/data directory

setwd(dir = "C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/Project3/manuscript/code/data")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attr <- listAttributes(ensembl)
rna_protein_matched <- read.csv("Fig2/RNA_protein_matched.csv")

gene_list <- unique(rna_protein_matched$gene_name)
res_ensembl_canonical <- getBM(attributes=c("external_gene_name","transcript_is_canonical","transcript_length","percentage_gene_gc_content","gene_biotype","ensembl_transcript_id","transcript_mane_select"), filters = "external_gene_name", values = gene_list, mart= ensembl)
# Only select canonical trnascripts:
res_ensembl_canonical <- res_ensembl_canonical %>%
  filter(transcript_is_canonical == 1)

res_ensembl_canonical <- res_ensembl_canonical %>%
  group_by(external_gene_name) %>%
  mutate(n_transcripts = n()) %>% # count transcripts per gene
  filter(n_transcripts == 1 | transcript_mane_select != "") %>%  # keep if only one row, or if it's the MANE Select transcript
  dplyr::select(-n_transcripts) %>% # remove helper column
  ungroup()

# Select transcripts: 
genes_transcripts <-res_ensembl_canonical$ensembl_transcript_id
res_ensembl_from_trscrpt <- getBM(attributes=c("external_gene_name","ensembl_transcript_id","transcript_length","cds_length"), filters = "ensembl_transcript_id", values = genes_transcripts, mart= ensembl)

#go from transcript to next page +
# use genomic features package

genome <- BSgenome.Hsapiens.UCSC.hg38
txdb <- makeTxDbFromGFF("Fig4/Homo_sapiens.GRCh38.114.gtf", format = "gtf")
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
utr3 <- threeUTRsByTranscript(txdb, use.names = TRUE)
cds <- cdsBy(txdb, use.names = TRUE)
all_tx <- transcriptsBy(txdb, by = "gene")

res_ensembl_from_trscrpt$utr3_seq <- NA
res_ensembl_from_trscrpt$utr5_seq <- NA
res_ensembl_from_trscrpt$cds_seq <- NA
res_ensembl_from_trscrpt$full_transcript_seq <- NA

for (i in 1:nrow(res_ensembl_from_trscrpt)) {
  trans <- res_ensembl_from_trscrpt$ensembl_transcript_id[i]
  gene <- res_ensembl_from_trscrpt$external_gene_name[i]
  # 5' UTR
  if (!is.null(utr5[[trans]])) {
    utr5_gr <- utr5[[trans]]
    seqlevelsStyle(utr5_gr) <- "UCSC"
    utr5_seq <- suppressWarnings(getSeq(genome, utr5_gr))
    full_utr5 <- paste(as.character(utr5_seq), collapse = "")
    res_ensembl_from_trscrpt$utr5_seq[i] <- full_utr5
  } else {
    res_ensembl_from_trscrpt$utr5_seq[i] <- NA
  }
  # 3' UTR
  if (!is.null(utr3[[trans]])) {
    utr3_gr <- utr3[[trans]]
    seqlevelsStyle(utr3_gr) <- "UCSC"
    utr3_seq <- suppressWarnings(getSeq(genome, utr3_gr))
    full_utr3 <- paste(as.character(utr3_seq), collapse = "")
    res_ensembl_from_trscrpt$utr3_seq[i] <- full_utr3
  } else {
    res_ensembl_from_trscrpt$utr3_seq[i] <- NA
  }
  # cds
  if (!is.null(cds[[trans]])) {
    cds_gr <- cds[[trans]]
    seqlevelsStyle(cds_gr) <- "UCSC"
    cds_seq <- suppressWarnings(getSeq(genome, cds_gr))
    full_cds <- paste(as.character(cds_seq), collapse = "")
    res_ensembl_from_trscrpt$cds_seq[i] <- full_cds
  } else {
    res_ensembl_from_trscrpt$cds_seq[i] <- NA
  }
  # full transcript
  if (!is.null(all_tx[[gene]])) {
    tx_gr <- all_tx[[gene]]
    seqlevelsStyle(tx_gr) <- "UCSC"
    tx_seq <- suppressWarnings(getSeq(genome, tx_gr))
    full_tx <- paste(as.character(tx_seq), collapse = "")
    res_ensembl_from_trscrpt$full_transcript_seq[i] <- full_tx
  } else {
    res_ensembl_from_trscrpt$full_transcript_seq[i] <- NA
  }
}

# Calculate GC content:
# Function to calculate GC content
gc_content <- function(seq) {
  if (is.na(seq) || nchar(seq) == 0) return(NA)
  seq <- toupper(seq)
  gc <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
  return((gc / nchar(seq)) * 100)
}

# Add columns
res_ensembl_from_trscrpt$utr5_length <- nchar(res_ensembl_from_trscrpt$utr5_seq)
res_ensembl_from_trscrpt$utr3_length <- nchar(res_ensembl_from_trscrpt$utr3_seq)
res_ensembl_from_trscrpt$cds_length_seq <- nchar(res_ensembl_from_trscrpt$cds_seq)

res_ensembl_from_trscrpt$utr5_gc <- sapply(res_ensembl_from_trscrpt$utr5_seq, gc_content)
res_ensembl_from_trscrpt$utr3_gc <- sapply(res_ensembl_from_trscrpt$utr3_seq, gc_content)
res_ensembl_from_trscrpt$cds_gc  <- sapply(res_ensembl_from_trscrpt$cds_seq, gc_content)

## Calculate the codon adaptation index (CAI):
# From: https://mt1022.github.io/cubar/articles/cubar.html
# Load transcriptomic data:
plot_dist <- function(x, xlab = 'values'){
  x <- stack(x)
  ggplot(x, aes(x = values)) +
    geom_histogram(bins = 40, fill = '#88CCEE') +
    labs(x = xlab, y = 'Number of genes') +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = 'black'))
}

rna <- read.csv(file = "transcriptomics/outputs/amol_per_ug_RNA.csv") # Remove all samples for prefiltered files (after qc eval)
# FIlter for matchign genes:
rna <- rna %>%
  filter(gene_name %in% res_ensembl_from_trscrpt$external_gene_name)
# Take median expression to select expression across all: 
rna$median <- apply(rna[, 1:32], 1, median)
heg <- head(rna[order(-rna$median), ], n = 200)
res_ensembl_from_trscrpt <- res_ensembl_from_trscrpt %>%filter(!is.na(cds_seq))
cf_all <- count_codons(res_ensembl_from_trscrpt$cds_seq)
rownames(cf_all) <- res_ensembl_from_trscrpt$external_gene_name

cf_all <- as.data.frame(cf_all) %>%
  tibble::rownames_to_column(var = "external_gene_name")
cf_heg <- cf_all %>%
  filter(external_gene_name %in% heg$gene_name)

rscu_heg <- est_rscu(cf_heg %>% dplyr::select(-external_gene_name))
rownames(cf_all) <- cf_all$external_gene_name
cf_all <- cf_all %>% dplyr::select(-external_gene_name) %>% as.matrix()
cai <- get_cai(cf_all, rscu_heg)
plot_dist(cai,'CAI')
df_cai <- data.frame(external_gene_name = names(cai),cai = cai)

# tAI calculation:
# From https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-summary.html
human_trna <- Biostrings::readRNAStringSet('Fig4/hg38-mature-tRNAs.fa')
# get tRNA gene copy number from GtRNADB

trna_gcn <- table(data.table::tstrsplit(sub(' .*', '', names(human_trna)), '-')[[3]])
trna_gcn <- trna_gcn[names(trna_gcn) != 'NNN'] # copy of each anticodon

ctab <- get_codon_table(gcid = '1')# Get codon table
trna_w <- est_trna_weight(trna_level = trna_gcn, codon_table = ctab)
tai <- get_tai(cf_all, trna_w = trna_w)
plot_dist(tai,'tAI')
hist(tai)
df_tai <- data.frame(external_gene_name = names(tai),tai = tai)

res_ensembl_from_trscrpt <- left_join(res_ensembl_from_trscrpt, df_cai, by = "external_gene_name")
res_ensembl_from_trscrpt <- left_join(res_ensembl_from_trscrpt, df_tai, by = "external_gene_name")

row.names(res_ensembl_from_trscrpt) <- res_ensembl_from_trscrpt$external_gene_name
res_ensembl_from_trscrpt$cds_seq[res_ensembl_from_trscrpt$external_gene_name == "17673"]

## Half-life analysis: half life from MCF7 cells
mrna_half_life <- read.csv("Fig4/GSE49831_MCF7_halflife.txt", sep = ' ')
colnames(mrna_half_life) <- c('MCF7_rep1','MCF7_rep2')
mrna_half_life$half_life <- apply(mrna_half_life[, 1:2], 1, mean)

hist(mrna_half_life$half_life,ylab = 'Number of genes',xlab = 'Half-life (in min)', breaks = 100)
mrna_half_life$external_gene_name <- row.names(mrna_half_life)
mrna_half_life <- mrna_half_life %>% dplyr::select(-MCF7_rep1, -MCF7_rep2)

# join tables to res_ensembl_from_trscrpt:
res_ensembl_from_trscrpt <- left_join(res_ensembl_from_trscrpt,mrna_half_life, by = "external_gene_name")

## Find nr of uORFs per gene:
fiveUTRs <- loadRegion(txdb, "leaders")

seqlevelsStyle(fiveUTRs) <- "UCSC"
seqlevels(fiveUTRs)[!seqlevels(fiveUTRs) %in% seqlevels(genome)] # not all chromosomes are in the genome assembly used
fiveUTRs <- keepSeqlevels(fiveUTRs, intersect(seqlevels(fiveUTRs), seqlevels(genome)), pruning.mode = "coarse")

tx_seqs <- extractTranscriptSeqs(genome, 
                                 fiveUTRs) 

fiveUTR_ORFs <- findMapORFs(fiveUTRs, tx_seqs, startCodon = startDefinition(1), stopCodon = stopDefinition(1),groupByTx = FALSE)
# fiveUTR_ORFs <- findMapORFs(fiveUTRs, tx_seqs,groupByTx = FALSE)

#txNames(fiveUTR_ORFs[1:2])
orfs_number <- as.data.frame(sort(table(txNames(fiveUTR_ORFs)), decreasing = TRUE))
colnames(orfs_number) <- c('ensembl_transcript_id','uORFs')

res_ensembl_from_trscrpt <- left_join(res_ensembl_from_trscrpt,orfs_number,by = 'ensembl_transcript_id')
res_ensembl_from_trscrpt$uORFs[is.na(res_ensembl_from_trscrpt$uORFs)] <- 0

# Get protein sequence:
for (i in 1:nrow(res_ensembl_from_trscrpt)){
  cds_seq <- DNAString(res_ensembl_from_trscrpt$cds_seq[i])
  res_ensembl_from_trscrpt$aa_seq[i] <- as.character(translate(cds_seq, if.fuzzy.codon = 'solve'))
}
res_ensembl_from_trscrpt$protein_length <- nchar(res_ensembl_from_trscrpt$aa_seq)
# molecular weight of protein:
res_ensembl_from_trscrpt$mw <- sapply(res_ensembl_from_trscrpt$aa_seq, mw)

for (i in 1:nrow(res_ensembl_from_trscrpt)){
  aa_seq <- res_ensembl_from_trscrpt$aa_seq[i]
  clean_seq <- gsub("\\*", "", aa_seq)
  molw <- mw(clean_seq)
  res_ensembl_from_trscrpt$protein_mw[i] <- molw
}

aa_freq <- alphabetFrequency(aa_seq, as.prob = TRUE)
aa_freq_df <- as.data.frame(aa_freq)

## Protein half-life:
protein_half_life <- read.csv("Fig2/MCF10A_half_life_protein_PMID29904729.txt",sep = '\t')
protein_half_life <- protein_half_life %>% # select only proteins iwth high R2:
  filter(R2 > 0.85) %>%
  dplyr::select(Gene,protein_half_life = half)
  
res_ensembl_from_trscrpt <- left_join(res_ensembl_from_trscrpt, protein_half_life, by = c("external_gene_name" = "Gene"))

#write.csv(res_ensembl_from_trscrpt,'5_processedData/20250714_gene_characteristics.csv',row.names = F)
res_ensembl_from_trscrpt <- read.csv('5_processedData/20250714_gene_characteristics.csv')

df <- res_ensembl_from_trscrpt

aa_freq_list <- lapply(df$aa_seq, function(seq) {
  aas <- strsplit(seq, split = "")[[1]]
  counts <- table(aas)
  freqs <- counts / length(aas)
  
  all_aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  freq_full <- setNames(rep(0, length(all_aas)), all_aas)
  freq_full[names(freqs)] <- freqs
  
  return(freq_full)
})

aa_freq_df <- do.call(rbind, aa_freq_list)

# Rename columns to add '_freq' 
colnames(aa_freq_df) <- paste0(colnames(aa_freq_df), "_freq")
aa_freq_df <- as.data.frame(aa_freq_df)
res_ensembl_from_trscrpt <- cbind(res_ensembl_from_trscrpt, aa_freq_df)
res_ensembl_from_trscrpt <- res_ensembl_from_trscrpt %>% select(!'*_freq')
#write.csv(res_ensembl_from_trscrpt,"gene_features.csv",row.names = F)

## General histograms of data

hist(res_ensembl_from_trscrpt$uORFs,breaks = 35)
hist(res_ensembl_from_trscrpt$utr5_length, breaks = 50)
hist(res_ensembl_from_trscrpt$utr3_length, breaks = 50)
hist(res_ensembl_from_trscrpt$utr5_gc, breaks = 50)
hist(res_ensembl_from_trscrpt$utr3_gc, breaks = 50)
hist(res_ensembl_from_trscrpt$half.1,breaks = 50)
