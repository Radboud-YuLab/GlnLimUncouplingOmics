## Script that plots Fig 2E-F
# written by: Cyriel Huijer
# last edited: 19 June 2026

library(tidyverse)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(ggthemes)
library(multcompView)
library(cowplot)

# Set working directory to /code/data directory
setwd("")

# clear environment
rm(list=ls())

# Read files
protein_rna_cell <- read_excel("cell_protein_RNA_overview.xlsx")
subset_gln <- read.csv(file = "Fig2/RNA_protein_matched.csv")

# Format protein_rna_cell: 
protein_rna_cell$sample_name <- paste0(protein_rna_cell$condition,"_",protein_rna_cell$replicate)

# Rerun:
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

head(subset_gln)
head(protein_rna_cell)

# select glutamine conditions, mean is already calculated so just the first reps are good
protein_rna_cell <- protein_rna_cell %>%
  filter(sample_label %in% c("GlnGlc_2","GlnGlc_3","GlnGlc_4","GlnGlc_5","GlnGlc_6")) %>%
  dplyr::select(condition,mean_ng_cell_gr,mean_ng_RNA_cell)

colnames(protein_rna_cell) <- c("condition","ng_protein_cell","ng_RNA_cell")

subset_gln <- subset_gln %>%
  left_join(protein_rna_cell, by = c("sample_name" = "condition"))

subset_gln$normalized_rna <- subset_gln$rna_amol_per_cell / subset_gln$ng_RNA_cell
subset_gln$normalized_protein <- subset_gln$protein_fmol_per_cell / subset_gln$ng_protein_cell

genes_pc2 <- read.csv("Fig2/pc2_genes_heatmap.csv")
genes_pc2 <- genes_pc2 %>%
  filter(gene %in% subset_gln$gene_name) # select only genes that are common in protein
genes_pc2$abs_pc2 <- abs(genes_pc2$PC2) 
genes_pc2 <- genes_pc2 %>% arrange(desc(abs_pc2))

subset_gln <- subset_gln %>% 
  filter(gene_name %in% genes_pc2$gene)
subset_gln$sample_name <- factor(subset_gln$sample_name, levels = c("Gln_04_mM", "Gln_12_mM","Gln_20_mM", "Gln_28_mM","Full_Med"))
plot_df <- subset_gln %>%
  dplyr::select(gene_name, sample_name, normalized_rna, normalized_protein) %>%
  pivot_longer(cols = c(normalized_rna, normalized_protein),
               names_to = "molecule", values_to = "value") %>%
  mutate(molecule = ifelse(molecule == "normalized_rna", "RNA", "Protein"))


# Dual axis pltos:
# make directory
#dir.create("gene_plots_dualaxis_scaled", showWarnings = FALSE)

for(g in unique(subset_gln$gene_name)) {
  
  df <- subset_gln %>% filter(gene_name == g)
  
  # scaling factor
  scale_factor <- max(df$normalized_rna, na.rm = TRUE) / 
    max(df$normalized_protein, na.rm = TRUE)
  
  # limits
  max_y <- max(c(df$normalized_rna, df$normalized_protein * scale_factor), na.rm = TRUE)
  
  p <- ggplot(df, aes(x = sample_name, group = 1)) +
    geom_line(aes(y = normalized_rna, color = "RNA"), size = 1.2) +
    geom_point(aes(y = normalized_rna, color = "RNA"), size = 3) +
    geom_line(aes(y = normalized_protein * scale_factor, color = "Protein"), size = 1.2) +
    geom_point(aes(y = normalized_protein * scale_factor, color = "Protein"), size = 3) +
    scale_y_continuous(
      name = "Normalized RNA",
      limits = c(0, max_y),
      sec.axis = sec_axis(~./scale_factor, name = "Normalized Protein")
    ) +
    scale_color_manual(values = c("RNA" = "lightgrey", "Protein" = "black")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          legend.position = "top",
          panel.grid.minor = element_blank()) +
    labs(title = g, x = "Condition")
  
  # ggsave(filename = paste0("gene_plots_dualaxis_scaled/", g, ".pdf"),
  #        plot = p, width = 4, height = 6)
}

# Pick YBX1 RPL35, ITGA2

# Individual plots:

label_map <- c(
  "Gln_04_mM" = "Gln (0.04 mM)",
  "Gln_12_mM" = "Gln (0.12 mM)",
  "Gln_20_mM" = "Gln (0.20 mM)",
  "Gln_28_mM" = "Gln (0.28 mM)",
  "Full_Med" = "Full Media (2.5 mM)"
)

plot_gene <- function(gene_to_plot) {
  df <- subset_gln %>% filter(gene_name == gene_to_plot)
  
  # Scaling factor for dual axis
  scale_factor <- max(df$normalized_rna, na.rm = TRUE) / 
    max(df$normalized_protein, na.rm = TRUE)
  
  max_y <- max(c(df$normalized_rna, df$normalized_protein * scale_factor), na.rm = TRUE)
  
  ggplot(df, aes(x = sample_name, group = 1)) +
    geom_line(aes(y = normalized_rna, color = "RNA"), size = 1.2) +
    geom_point(aes(y = normalized_rna, color = "RNA"), size = 3) +
    geom_line(aes(y = normalized_protein * scale_factor, color = "Protein"), size = 1.2) +
    geom_point(aes(y = normalized_protein * scale_factor, color = "Protein"), size = 3) +
    scale_y_continuous(
      name = "Normalized RNA",
      limits = c(0, max_y),
      sec.axis = sec_axis(~./scale_factor, name = "Normalized Protein")
    ) +
    scale_x_discrete(labels = label_map) +   # <-- apply mapping here
    scale_color_manual(values = c("RNA" = "lightgrey", "Protein" = "black")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(),
          legend.position = "top") +
    labs(title = gene_to_plot, x = "")
}


# Generate plots for ITGA2, ENO1, FASN (Fig2E-G)
p1 <- plot_gene("YBX1")
p2 <- plot_gene("RPL35")
p3 <- plot_gene("ITGA2")

# Arrange side by side or in grid
ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "top")

