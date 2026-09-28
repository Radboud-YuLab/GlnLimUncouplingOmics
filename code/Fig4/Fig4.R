## Gene features script (Fig 4)
# Written by: Cyriel Huijer
# Last edited: 19 June 2026

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

# Load in files

gene_corr_gln <- read.csv(file = "Fig4/gene_correlations.csv")
gene_data <- read.csv(file = "gene_features.csv")

gene_corr_gln
shapiro.test(gene_corr_gln$pearson_corr)

#they are  are non normally distributed, so lowest % genes and top % genes are picked for further analysis

# Join gene_corr_gln with gene_data
gene_corr_gln <- gene_corr_gln %>%
  left_join(gene_data, by = c("gene_name" = "external_gene_name"))

gene_corr_gln_backup <- gene_corr_gln

gene_corr_gln <- gene_corr_gln %>% filter(corr_group %in% c("low","high"))
gene_corr_gln$corr_group <- factor(gene_corr_gln$corr_group, levels = c("low", "high"))


# Bubble plot:

numeric_vars <- c("transcript_length", "cds_length", "utr5_length", "utr3_length",
                  "utr5_gc", "utr3_gc", "cds_gc", "cai", "tai", "half_life",
                  "mw", "protein_half_life",
                  "A_freq", "C_freq", "D_freq", "E_freq", "F_freq", "G_freq", "H_freq", "I_freq", 
                  "K_freq", "L_freq", "M_freq", "N_freq", "P_freq", "Q_freq", "R_freq", "S_freq",
                  "T_freq", "V_freq", "W_freq", "Y_freq")

# Wilcoxon tests and calculate log2FC
results <- lapply(numeric_vars, function(var) {
  wt <- wilcox.test(gene_corr_gln[[var]] ~ gene_corr_gln$corr_group)
  
  med_high <- median(gene_corr_gln[[var]][gene_corr_gln$corr_group == "high"], na.rm = TRUE)
  med_low  <- median(gene_corr_gln[[var]][gene_corr_gln$corr_group == "low"], na.rm = TRUE)
  
  log2FC   <- log2(med_high / med_low)
  perc_diff <- ((med_high - med_low) / med_low) * 100   # <- percentage difference
  
  tibble(
    variable = var,
    med_high = med_high,
    med_low = med_low,
    log2FC   = log2FC,
    perc_diff = perc_diff,
    p.value  = wt$p.value
  )
}) %>% bind_rows() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         significance = -log10(p.adj),
         sig_label = p.adj < 0.05)

# Fig 4A-C plots
results_4A <- results %>% filter(variable %in% c("utr3_length","transcript_length","cds_length","mw","utr5_length","utr3_gc","cds_gc","utr5_gc"))
ggplot(results_4A, aes(x = perc_diff, y = reorder(variable, perc_diff))) +
  geom_point(size = 5,aes(
    color = significance,
    shape = sig_label
  )) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#006CD1", high = "#994F00") +
  labs(
    x = "Percentage difference (high vs low correlating genes)",
    y = "Feature",
    color = "-log10(adj. p)",
    shape = "Significant"
  ) +
  theme_bw(base_size = 14)
#ggsave("6_results/final/fig4/20260129_Fig4A.pdf", width = 6, height = 4)

results_4B <- results %>% filter(variable %in% c("half_life","cai","tai","protein_half_life"))
ggplot(results_4B, aes(x = perc_diff, y = reorder(variable, perc_diff))) +
  geom_point(size = 5,aes(
    color = significance,
    shape = sig_label
  )) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#006CD1", high = "#994F00") +
  labs(
    x = "Percentage difference (high vs low correlating genes)",
    y = "Feature",
    color = "-log10(FDR)",
    shape = "Significant"
  ) +
  theme_bw(base_size = 14)
#ggsave("6_results/final/fig4/20260129_Fig4B.pdf", width = 6, height = 4)

results_4C <- results %>% filter(grepl("freq$", variable))
ggplot(results_4C, aes(x = perc_diff, y = reorder(variable, perc_diff))) +
  geom_point(size = 5,aes(
    color = significance,
    shape = sig_label
  )) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "#006CD1", high = "#994F00") +
  labs(
    x = "Percentage difference (high vs low correlating genes)",
    y = "Feature",
    color = "-log10(adj. p)",
    shape = "Significant"
  ) +
  theme_bw(base_size = 14)
#ggsave("6_results/final/fig4/20260129_Fig4C.pdf", width = 6, height = 8.64)


# Individual plots of features
plot_feature_boxplot <- function(data, feature, condition = "Gln") {
  # convert feature name to symbol
  feat_sym <- rlang::sym(feature)
  # run t-test
  wilcox_res <- wilcox.test(data[[feature]] ~ data$corr_group)
  print(wilcox_res)
  # make plot
  p <- ggplot(data, aes(x = corr_group, y = !!feat_sym, fill = corr_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    labs(
      x = "Correlation group",
      y = feature,
      title = paste(feature, " in high vs low correlating genes (", condition, " limitation)", sep = "")
    ) +
    scale_fill_manual(values = c("low" = "tomato", "high" = "steelblue")) +
    theme_bw(base_size = 14) +
    stat_compare_means(method = "wilcox.test", label = "p.format")
  
  return(p)
}

# Individual plots
# plot_feature_boxplot(gene_corr_gln, "utr5_length", "Gln")
# plot_feature_boxplot(gene_corr_gln, "cai", "Gln")
# plot_feature_boxplot(gene_corr_gln, "tai", "Gln")
# plot_feature_boxplot(gene_corr_gln, "cds_length", "Gln")
# plot_feature_boxplot(gene_corr_gln, "utr3_length", "Gln")
# plot_feature_boxplot(gene_corr_gln, "utr5_gc", "Gln")
# plot_feature_boxplot(gene_corr_gln, "utr3_gc", "Gln")
# plot_feature_boxplot(gene_corr_gln, "cds_gc", "Gln")
# plot_feature_boxplot(gene_corr_gln, "transcript_length", "Gln")
# plot_feature_boxplot(gene_corr_gln, "half_life", "Gln")
# plot_feature_boxplot(gene_corr_gln, "protein_half_life", "Gln")
# plot_feature_boxplot(gene_corr_gln, "uORFs", "Gln")
# plot_feature_boxplot(gene_corr_gln, "protein_length", "Gln")
