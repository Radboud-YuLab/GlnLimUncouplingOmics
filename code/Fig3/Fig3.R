## This script takes mol/ug RNA/protein to mol/cell for all conditions:
# Written by: Cyriel Huijer
# Date: 27 May 2025
# Last edit: 19 June 2026

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
setwd("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/Project3/manuscript/code/data")

# clear environment
rm(list=ls())

## Format RNA/protein data:

# Load in abs quant data:
protein <- read.csv(file = "protein_fmol_cell.csv")
row.names(protein) <- protein$gene_name
protein <- protein %>% dplyr::select(!gene_name)
rna <- read.csv(file = "rna_amol_cell.csv")
rna <- rna %>% filter(!is.na(gene_name))
row.names(rna) <- rna$gene_name
rna <- rna %>% dplyr::select(!gene_name)

# select conditions after manual QC
rna <- rna %>% dplyr::select(-c(Glc_00_gL_1,Glc_70_gL_1,Glc_70_gL_2,Glc_70_gL_3,Gln_00_mM_3,Gln_12_mM_1,Gln_28_mM_3,Gln_20_mM_3,Gln_00_mM_2,Gln_00_mM_1)) # PCA is good for Gln
protein <- protein %>% dplyr::select(!ends_with("2")) %>% dplyr::select(-c(Gln_12_mM_1))

conditions <- rna %>% 
  colnames() %>% 
  str_remove("_\\d+$") %>% 
  unique()
# Calculate mean per condition across replicates
rna_avg <- sapply(conditions, function(cond) {
  cols <- grep(paste0("^", cond, "_\\d+$"), colnames(rna), value = TRUE)
  rowMeans(rna[, cols, drop = FALSE])
}) %>% 
  as.data.frame()
# Calculate mean expression across replicates per condition
protein_avg <- sapply(conditions, function(cond) {
  cols <- grep(paste0("^", cond, "_\\d+$"), colnames(protein), value = TRUE)
  rowMeans(protein[, cols, drop = FALSE])
}) %>%
  as.data.frame()

protein_avg <- protein_avg[apply(protein_avg, 1, function(x) all(x > 1e-8)), ]


rna_long <- rna_avg %>%
  mutate(gene_name = rownames(.)) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample_name",
    values_to = "rna_amol_per_cell"
  )

protein_long <- protein_avg %>%
  mutate(gene_name = rownames(.)) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample_name",
    values_to = "protein_fmol_per_cell"
  )

protein_long <- protein_long %>%
  filter(protein_fmol_per_cell > 1e-8)

# Match samples:

matched_samples <- intersect(rna_long$sample_name, protein_long$sample_name)

rna_long <- rna_long %>%
  filter(sample_name %in% matched_samples) %>%
  dplyr::select(gene_name, sample_name, rna_amol_per_cell)
#write_csv(rna_long,"Fig2/all_RNA_heatmap.csv")

protein_long <- protein_long %>%
  filter(sample_name %in% matched_samples) %>%
  dplyr::select(gene_name, sample_name, protein_fmol_per_cell)

# Join by gene_name and sample_name (only exact matches)
rna_protein_matched <- inner_join(rna_long, protein_long,
                                  by = c("gene_name", "sample_name"))

#write_csv(rna_protein_matched,"5_processedData/20250917_RNA_protein_matched.csv")

# Make subsets
subset_gln <- rna_protein_matched %>%
  filter(str_detect(sample_name, "Gln") | str_detect(sample_name, "Full_Med"))
# for fig 2E-G
#write.csv(subset_gln,"Fig2/RNA_protein_matched.csv",row.names = F)
rna_protein_matched <- rna_protein_matched %>% filter(str_detect(sample_name, "Gln") | str_detect(sample_name, "Full_Med") | str_detect(sample_name, "Glc_06_gL|Glc_29_gL|Glc_18_gL"))

## Figure 3A: --------------------------------------------------------------------

rna_protein_matched <- rna_protein_matched %>%
  filter(sample_name == "Full_Med") %>%
  filter(log10(protein_fmol_per_cell) > -6)
model <- lm(log10(rna_protein_matched$protein_fmol_per_cell) ~ log10(rna_protein_matched$rna_amol_per_cell))
summary(model)
r_val <- cor(log10(rna_protein_matched$rna_amol_per_cell),
             log10(rna_protein_matched$protein_fmol_per_cell),
             method = "pearson")

ggplot(rna_protein_matched, aes(x = rna_amol_per_cell, y = protein_fmol_per_cell)) +
  geom_point(color = "black") +
  theme_bw(base_size = 15) +
  annotate(
    "text",
    x = min(rna_protein_matched$rna_amol_per_cell, na.rm = TRUE),
    y = max(rna_protein_matched$protein_fmol_per_cell, na.rm = TRUE),
    label = paste0("r = ", round(r_val, 2)),
    hjust = 0, vjust = 1.5, size = 5
  ) +
  ylab("Protein abundance (proteins/cell)") +
  xlab("RNA abundance (transcripts/cell)") +
  scale_x_log10(labels = scales::label_log()) +
  scale_y_log10(labels = scales::label_log())
#ggsave(filename = "6_results/final/fig3/dynamic_range_full_med_molecules.pdf", width = 5, height = 5)

## Fig 3B: Within-gene correlation ---------------------------------------------

# Function to calculate Pearson correlation per gene
calculate_gene_correlation <- function(df) {
  df %>%
    group_by(gene_name) %>%
    summarize(
      n = n(),
      spearman_corr = cor(rna_amol_per_cell, protein_fmol_per_cell, method = "spearman"),
      pearson_corr = cor(log(rna_amol_per_cell), log(protein_fmol_per_cell), method = "pearson"),
      .groups = "drop"
    ) %>%
    filter(n >= 2)  # Only keep genes with at least 2 paired samples
}

# Run correlation per gene
gene_corr_gln <- calculate_gene_correlation(subset_gln)
gene_corr_all <- calculate_gene_correlation(rna_protein_matched)

#Add labels for Gln
q_gln <- quantile(gene_corr_gln$pearson_corr, probs = c(0.25, 0.75), na.rm = TRUE)

# Figure 3B
ggplot(gene_corr_gln, aes(x = pearson_corr)) +
  geom_histogram(aes(y = ..count..),binwidth = 0.1, color = "black",alpha = 0.5) +  
  geom_vline(xintercept = median(gene_corr_gln$pearson_corr, na.rm = TRUE), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = q_gln[1], linetype = "dashed", color = "black",size = 1) +# 25 percentile
  geom_vline(xintercept = q_gln[2], linetype = "dashed", color = "black",size = 1) +  # 75 percentile
  #geom_density()+
  labs(
    title = "",
    x = "Pearson Correlation",
    y = "Count"
  ) +
  theme_bw(base_size = 15)+
  xlim(-1,1)
#ggsave("6_results/final/fig3/histogram_Gln.pdf", width=5, height=5)

median(gene_corr_gln$pearson_corr,na.rm = T) # 0.36
#median(gene_corr_gln_filter$pearson_corr,na.rm = T)

## Fig 3C-D: ------------------------------------------

#Add labels for Gln
q_gln <- quantile(gene_corr_gln$pearson_corr, probs = c(0.25, 0.75), na.rm = TRUE)
gene_corr_gln <- gene_corr_gln %>%
  mutate(corr_group = case_when(
    pearson_corr <= q_gln[1] ~ "low",
    pearson_corr >= q_gln[2] ~ "high",
    TRUE ~ "mid"
  ))

# Write correlation dataframe for Fig4
write.csv(gene_corr_gln,"Fig4/gene_correlations.csv")
# Glutamine

# 1. Filter for correlated genes
correlated_genes <- gene_corr_gln %>%
  filter(corr_group == "low") %>%
  pull(gene_name)

# 2. Run GO enrichment
GO_corr_results_low <- enrichGO(
  gene         = correlated_genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Convert results to dataframe
GO_corr_df_low <- as.data.frame(GO_corr_results_low)
GO_corr_df_low$quartile <- "low"

# 1. Filter for correlated genes
correlated_genes <- gene_corr_gln %>%
  filter(corr_group == "high") %>%
  pull(gene_name)

# 2. Run GO enrichment
GO_corr_results_high <- enrichGO(
  gene         = correlated_genes,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Convert results to dataframe
GO_corr_df_high <- as.data.frame(GO_corr_results_high)
GO_corr_df_high$quartile <- "high"


GO_corr_df <- rbind(GO_corr_df_high,GO_corr_df_low)
#write.csv(GO_corr_df, "../../supplementary_material/TableS1_2.csv",row.names = F)

# 3. Plot top GO terms

# high correlation genes
# Take top 10 GO terms
top_GO_corr <- GO_corr_df_high %>%
  slice_head(n = 10) %>%
  mutate(
    # Convert GeneRatio to fraction
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    # Rank order by GeneRatio (largest at top)
    Description = factor(Description, levels = Description[order(GeneRatio_num, decreasing = FALSE)])
  )

# Bubble plot
p1 <-ggplot(top_GO_corr, aes(x = Description, y = GeneRatio_num)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count") +
  coord_flip() +
  labs(
    x = "",
    y = "Gene Ratio",
    title = "Low RNA–Protein Correlating Genes",
    subtitle = "Biological Process (BP)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text()
  )

p1

# Low correlation genes
# Take top 10 GO terms
top_GO_corr <- GO_corr_df_low %>%
  slice_head(n = 10) %>%
  mutate(
    # Convert GeneRatio to fraction
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    # Rank order by GeneRatio (largest at top)
    Description = factor(Description, levels = Description[order(GeneRatio_num, decreasing = FALSE)])
  )

# Bubble plot
p2 <-ggplot(top_GO_corr, aes(x = Description, y = GeneRatio_num)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count") +
  coord_flip() +
  labs(
    x = "",
    y = "Gene Ratio",
    title = "Low RNA–Protein Correlating Genes",
    subtitle = "Biological Process (BP)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text()
  )

p2
aligned <- align_plots(p1, p2, align = "v", axis = "l")


# Figure 3C-D
plot_grid(aligned[[1]], aligned[[2]], ncol = 1)

