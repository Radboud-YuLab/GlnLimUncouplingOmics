## This script takes mol/ug RNA/protein to mol/cell for all conditions:
# Written by: Cyriel Huijer
# Date: 27 May 2025

library(tidyverse)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(ggthemes)
library(multcompView)
library(cowplot)

# clear environment
rm(list=ls())

# Set working directory to /code/data directory
setwd(dir = "")

# Load in data:
rna <- read.csv(file = "transcriptomics/outputs/amol_per_ug_RNA.csv") 
protein <- read.csv(file = "proteomics/output/fmol_ug_protein.csv")
colnames(protein) <- gsub("_fmol", "", colnames(protein))
protein_rna_cell <- read_excel("cell_protein_RNA_overview.xlsx")

# Format protein_rna_cell: 
protein_rna_cell$sample_name <- paste0(protein_rna_cell$condition,"_",protein_rna_cell$replicate)

# Format:
protein_rna_cell <- protein_rna_cell %>%
  left_join(
    protein_rna_cell %>%
      group_by(condition) %>%
      summarise(
        mean_ng_cell_gr = mean(ng_cell_gr[replicate %in% c(1, 3)], na.rm = TRUE),
        mean_ng_RNA_cell = mean(ng_RNA_cell, na.rm = TRUE)
      ),
    by = "condition"
  )

# Determine overlap in RNA/protein data:
length(intersect(rna$gene_name,protein$gene_name)) # 1705 genes are overlapping between both datasets

# Calculate protein/cell and rna/cell:

# Pivot longer: each row = gene-sample pair
rna_long <- rna %>%
  pivot_longer(cols = -gene_name,names_to = "sample_name",values_to = "rna_amol_ug")

protein_long <- protein %>%
  pivot_longer(cols = -gene_name, names_to = "sample_name", values_to = "protein_fmol_ug")

# Merge dataframes:

rna_long <- rna_long %>%
  left_join(protein_rna_cell %>% dplyr::select(sample_name, ng_RNA_cell, mean_ng_RNA_cell), by = "sample_name")

protein_long <- protein_long %>%
  left_join(protein_rna_cell %>% dplyr::select(sample_name, ng_cell_gr, mean_ng_cell_gr), by = "sample_name")

# convert to rna/cell or protein/cell:

rna_long <- rna_long %>%
  mutate(rna_amol_per_cell = rna_amol_ug * (mean_ng_RNA_cell / 1000)) # this is for the mean RNA per condition
protein_long <- protein_long %>%
  mutate(protein_fmol_per_cell = protein_fmol_ug * (mean_ng_cell_gr / 1000)) # this is for the mean protien per condition

# For export:

rna_long <- rna_long %>%
  dplyr::select(gene_name, sample_name, rna_amol_per_cell)
# Remove duplicates (keep first)
rna_long_unique <- rna_long %>%
  distinct(gene_name, sample_name, .keep_all = TRUE)
# Pivot to wide format
rna_wide <- rna_long_unique %>%
  pivot_wider(names_from = sample_name, values_from = rna_amol_per_cell)
rna_wide <- rna_wide %>% filter(!is.na(gene_name))
rna_wide <- rna_wide %>% filter(rowSums(dplyr::select(., -gene_name)) != 0)

#write.csv(rna_wide,file = "rna_amol_cell.csv",row.names = F)

protein_long <- protein_long %>%
  dplyr::select(gene_name, sample_name, protein_fmol_per_cell)
# Remove duplicates (keep first)
protein_long_unique <- protein_long %>%
  distinct(gene_name, sample_name, .keep_all = TRUE)
# Pivot to wide format
protein_wide <- protein_long_unique %>%
  pivot_wider(names_from = sample_name, values_from = protein_fmol_per_cell)
#write.csv(protein_wide,file = "protein_fmol_cell.csv",row.names = F)

