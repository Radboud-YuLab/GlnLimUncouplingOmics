## Script for correlation PCA1 with transcript/protein per cell plots for transcriptomics and proteomics on absquant data (mol/cell)
# Written by: Cyriel Huijer
# Date: 03 June 2025
# Last updated: 19 June 2026

library(ggrepel)
library(tidyverse)
library(readxl)

# Set working directory:
setwd("")

# clear environment
rm(list=ls())

# Load in abs quant data:
protein <- read.csv(file = "data/protein_fmol_cell.csv")
row.names(protein) <- protein$gene_name
protein <- protein %>% dplyr::select(!gene_name)
rna <- read.csv(file = "data/rna_amol_cell.csv")
rna <- rna %>% filter(!is.na(gene_name))
row.names(rna) <- rna$gene_name
rna <- rna %>% dplyr::select(!gene_name)

length(intersect(row.names(protein),row.names(rna))) # 1705 genes are overlapping between both datasets

# RNAseq ----------------------------------------------------------------------
# Select glutamine conditions: 

gln_rna <- rna %>% dplyr::select(contains("Gln"),contains("Full_Med")) %>%filter(apply(., 1, var) != 0)
gln_rna <- gln_rna %>% dplyr::select(-c(Gln_00_mM_3,Gln_12_mM_1,Gln_28_mM_3,Gln_20_mM_3,Gln_00_mM_2,Gln_00_mM_1)) %>%filter(apply(., 1, var) != 0)

gln_long <- gln_rna %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(condition = sub("_(\\d+)$", "", sample))

# Average across replicates
gln_avg <- gln_long %>%
  group_by(gene, condition) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expr) %>%
  column_to_rownames("gene")

gln_avg_rna <- gln_avg

## Get PCA PC2 highest loadings for Figure 2B (heatmap)
# Run PCA on averaged expression
pca_res <- prcomp(t(gln_avg), scale. = TRUE)

# Extract loadings (rotation tells how much each gene contributes to each PC)
loadings <- as.data.frame(pca_res$rotation)

# Add gene names
loadings$gene <- rownames(loadings)

# Order by absolute contribution to PC2
top_pc2 <- loadings %>%
  filter(abs(PC2) >= 0.015) %>%
  arrange(desc(abs(PC2)))
#write_csv(x = top_pc2,"data/Fig2/pc2_genes_heatmap.csv")

# PC1 correlates with RNA/cell (Figure 2A)

custom_colors <- c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "#5F8247")
custom_labels <- c("Gln (0.04 mM)", "Gln (0.12 mM)", "Gln (0.20 mM)", "Gln (0.28 mM)", "Full Media (2.5 mM)")

protein_rna_cell <- read_excel("data/cell_protein_RNA_overview.xlsx")
protein_rna_cell$sample_name <- paste0(protein_rna_cell$condition,"_",as.character(protein_rna_cell$replicate))
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

protein_rna_cell <- protein_rna_cell %>%
  filter(sample_label %in% c("GlnGlc_2","GlnGlc_3","GlnGlc_4","GlnGlc_5","GlnGlc_6")) %>%
  dplyr::select(condition,mean_ng_cell_gr,mean_ng_RNA_cell)

gln_rna <- rna %>% dplyr::select(contains("Gln"),contains("Full_Med")) %>%filter(apply(., 1, var) != 0)
gln_rna <- gln_rna %>% dplyr::select(-c(Gln_00_mM_3,Gln_12_mM_1,Gln_28_mM_3,Gln_20_mM_3,Gln_00_mM_2,Gln_00_mM_1)) %>%filter(apply(., 1, var) != 0)
gln_long <- gln_rna %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(condition = sub("_(\\d+)$", "", sample))

# Average across replicates
gln_avg <- gln_long %>%
  group_by(gene, condition) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expr) %>%
  column_to_rownames("gene")

sample_data <- gln_avg %>% t() %>% as.data.frame()
colnames(sample_data) <- row.names(gln_avg)
sample_data$sample_name <- rownames(sample_data)
# Extract condition and replicate
# sample_data <- sample_data %>%
#   mutate(condition = sub("_(\\d+)$", "", sample_name), replicate = sub(".*_(\\d+)$", "\\1", sample_name))
sample_data <- sample_data %>%
  mutate(condition = sub("_(\\d+)$", "", sample_name))
# Perform PCA
pca_res <- prcomp(sample_data %>% dplyr::select(-sample_name, -condition), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>%
  mutate(condition = sample_data$condition)
pca_df$sample_name <- row.names(pca_df)

comb_table <- left_join(protein_rna_cell, pca_df, by = "condition")
comb_table <- comb_table[!is.na(comb_table$PC1),]
head(comb_table)

# Make sure your factor levels match the order of colors/labels
comb_table$condition <- factor(comb_table$condition,
                                    levels = c("Gln_04_mM", "Gln_12_mM", 
                                               "Gln_20_mM", "Gln_28_mM", 
                                               "Full_Med"))
label_map <- c(
  "Gln_04_mM" = "Gln (0.04 mM)",
  "Gln_12_mM" = "Gln (0.12 mM)",
  "Gln_20_mM" = "Gln (0.20 mM)",
  "Gln_28_mM" = "Gln (0.28 mM)",
  "Full_Med" = "Full Media (2.5 mM)"
)

comb_table$custom_label <- label_map[as.character(comb_table$condition)]

# Get % explained for plot: 
summ <- summary(pca_res)
pc1_pct <- round(summ$importance[2, 1] * 100, 1)  # multiply by 100 for %

ggplot(comb_table, aes(x = PC1, y = mean_ng_RNA_cell, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = "ng RNA/cell",
    color = "Condition",
    title = ""
  ) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") +
  geom_text(
    data = comb_table %>% distinct(condition, .keep_all = TRUE),
    aes(label = custom_label),
    vjust = -1, size = 3
  )
## Figure 2A:
# ggsave(
#   paste0("./6_results/final/fig2/", format(Sys.Date(), "%Y%m%d"), "_rna_PC1.pdf"),
#   width = 5, height = 5
# )
# Proteomics ------------------------------------------------------------------

# Select gLutamine samples

gln_protein <- protein %>% dplyr::select(contains("Gln"),contains("Full_Med")) %>%filter(apply(., 1, var) != 0)
gln_protein <- gln_protein %>% dplyr::select(!ends_with("2"),-c(Gln_00_mM_3,Gln_12_mM_1)) %>% filter(apply(., 1, var) != 0)

gln_long <- gln_protein %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  mutate(condition = sub("_(\\d+)$", "", sample))

# Average across replicates
gln_avg <- gln_long %>%
  group_by(gene, condition) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expr) %>%
  column_to_rownames("gene")
gln_avg_protein <- gln_avg

## DO the same for proteomics data: -------------------------------------------

sample_data <- gln_avg_protein %>% t() %>% as.data.frame()
colnames(sample_data) <- row.names(gln_protein)
sample_data$sample_name <- rownames(sample_data)
# Extract condition and replicate
sample_data <- sample_data %>%
  mutate(condition = sub("_(\\d+)$", "", sample_name))
# Perform PCA
pca_res <- prcomp(sample_data %>% dplyr::select(-sample_name, -condition), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x) %>%
  mutate(condition = sample_data$condition)
pca_df$sample_name <- row.names(pca_df)

comb_table <- left_join(protein_rna_cell, pca_df, by = "condition")
comb_table <- comb_table[!is.na(comb_table$PC1),]
head(comb_table)

# Make sure your factor levels match the order of colors/labels
comb_table$condition <- factor(comb_table$condition,
                               levels = c("Gln_04_mM", "Gln_12_mM", 
                                          "Gln_20_mM", "Gln_28_mM", 
                                          "Full_Med"))
comb_table$custom_label <- label_map[as.character(comb_table$condition)]

# Get % explained for plot: 
summ <- summary(pca_res)
pc1_pct <- round(summ$importance[2, 1] * 100, 1)  # multiply by 100 for %

ggplot(comb_table, aes(x = PC1, y = mean_ng_cell_gr , color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = F, color = "black", linetype = "dashed") +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = "ng protein/cell",
    color = "Condition",
    title = ""
  ) +
  scale_color_manual(values = custom_colors, labels = custom_labels) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") +
  geom_text(
    data = comb_table %>% distinct(condition, .keep_all = TRUE),
    aes(label = custom_label),
    vjust = -1, size = 3
  )

# # Fig 2C
# ggsave(
#   paste0("./6_results/final/fig2/", format(Sys.Date(), "%Y%m%d"), "_protein_PC1.pdf"),
#   width = 5, height = 5
# )

