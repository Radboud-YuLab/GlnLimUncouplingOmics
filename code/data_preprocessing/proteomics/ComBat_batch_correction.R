## SVA ComBat batch correction
# Written by: Cyriel Huijer
# Date: 29 April 2025

# Set working directory to /code/data directory
setwd("")

# Load in libraries:
library(tidyverse) # tidyverse_2.0.0
library(sva) # sva_3.46.0

# Clear workspace:

rm(list=ls())

# Use colorblind-friendly palette
cbPalette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999", "#000000", "#AAAA00", "#44AA99"
)


# Load in data:

norm_1 <- read.csv('proteomics/output/rep1_vsn_normalized_correctlabels.csv', sep = ',')
norm_2 <- read.csv('proteomics/output/rep2_vsn_normalized_correctlabels.csv', sep = ',')
norm_3 <- read.csv('proteomics/output/rep3_vsn_normalized_correctlabels.csv', sep = ',')

# Normalize for the pool

data_corr_1_ratio <- sweep(norm_1, 1, norm_1[, 1], FUN = "-")
data_corr_2_ratio <- sweep(norm_2, 1, norm_2[, 1], FUN = "-")
data_corr_3_ratio <- sweep(norm_3, 1, norm_3[, 2], FUN = "-")  # Use column 2 for the pooled sample for this dataframe

# Find common rownames
common_genes_ratio <- Reduce(intersect, list(rownames(data_corr_1_ratio),
                                             rownames(data_corr_2_ratio),
                                             rownames(data_corr_3_ratio)))
# Subset each dataset to only common genes
data_corr_1_ratio <- data_corr_1_ratio[common_genes_ratio, ]
data_corr_2_ratio <- data_corr_2_ratio[common_genes_ratio, ]
data_corr_3_ratio <- data_corr_3_ratio[common_genes_ratio, ]

combined_data_ratio <- cbind(data_corr_1_ratio, data_corr_2_ratio, data_corr_3_ratio)
pca_input <- t(combined_data_ratio)

pca_result <- prcomp(pca_input, scale. = TRUE)

# Make a data frame with PC scores and sample labels
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

# Still a clear batch effect:
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text_repel(size = 3, max.overlaps = 50) +
  theme_minimal() +
  ggtitle("") +
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "% variance)"))

## Do the batch correction with ComBat from sva package:

# convert combined data to matrix
expr_matrix <- as.matrix(combined_data_ratio)

# create batch vector dynamically from column names
batch_vector <- as.numeric(sub(".*_([123])$", "\\1", colnames(expr_matrix)))

# create condition vector dynamically from column names
condition_vector <- gsub("_[123]$", "", colnames(expr_matrix))

# remove all pooled samples and the entire second batch
exclude_samples <- !(condition_vector == "pooled" | batch_vector == 2 | (condition_vector == "Gln_12_mM" & batch_vector == 1))
#exclude_samples <- !(condition_vector == "pooled")


# apply filtering
expr_matrix_filtered <- expr_matrix[, exclude_samples]
batch_vector_filtered <- batch_vector[exclude_samples]
condition_vector_filtered <- condition_vector[exclude_samples]

# build model matrix for the remaining conditions
mod_filtered <- model.matrix(~ as.factor(condition_vector_filtered))

# apply ComBat batch correction
expr_matrix_corrected <- ComBat(dat = expr_matrix_filtered, batch = batch_vector_filtered, mod = mod_filtered)

# save expression matrix for further analysis:
expr_matrix_corrected_df <- data.frame(Accession = rownames(expr_matrix_corrected), expr_matrix_corrected)
#write.csv(expr_matrix_corrected_df, "proteomics/output/vsn_combat_corrected_log2_pool_ratios_correctlabels.csv", row.names = FALSE)

# perform PCA
pca_result <- prcomp(t(expr_matrix_corrected), scale. = TRUE)

# prepare PCA results for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4],
  Batch = factor(batch_vector_filtered),
  Condition = factor(condition_vector_filtered),
  Sample = colnames(expr_matrix_filtered)  # Add sample names
)

# plot PCA
explained_variance <- summary(pca_result)$importance[2, ] * 100  # Percent of variance explained

# PC1 vs PC2

ggplot(pca_df, aes(x = PC1, y = PC2, shape = Batch, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 50) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = cbPalette)+
  labs(
    title = "PCA After ComBat Correction (Batch 2 Removed)",
    x = paste("PC1 -", round(explained_variance[1], 2), "% Explained Variance"),
    y = paste("PC2 -", round(explained_variance[2], 2), "% Explained Variance")
  ) +
  theme(legend.position = "right")

# PC1 vs PC3

ggplot(pca_df, aes(x = PC1, y = PC3, shape = Batch, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 50) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = cbPalette)+
  labs(
    title = "PCA After ComBat Correction (Batch 2 Removed)",
    x = paste("PC1 -", round(explained_variance[1], 2), "% Explained Variance"),
    y = paste("PC3 -", round(explained_variance[3], 2), "% Explained Variance")
  ) +
  theme(legend.position = "right")

# PC1 vs PC4
ggplot(pca_df, aes(x = PC1, y = PC4, shape = Batch, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 50) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = cbPalette)+
  labs(
    title = "PCA After ComBat Correction (Batch 2 Removed)",
    x = paste("PC1 -", round(explained_variance[1], 2), "% Explained Variance"),
    y = paste("PC4 -", round(explained_variance[4], 2), "% Explained Variance")
  ) +
  theme(legend.position = "right")


# Select Gln only:
# Gln conditions:
gln_conditions <- expr_matrix_corrected_df %>%
  dplyr::select(matches("Gln_00|Gln_04|Gln_12|Gln_20|Full_Med"))
# PCA
pca_res <- prcomp(t(gln_conditions))

pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- colnames(gln_conditions)
pca_df$condition <- sub("_[^_]+$", "", pca_df$sample)

# Percent variance
percentVar <- round(100 * summary(pca_res)$importance[2, 1:5], 1)

# Plot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3,aes(color = condition)) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = cbPalette) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )+
  ggtitle("PCA plot of glutamine limitation samples (TMT proteomics)")
#ggsave("6_results/Proteomics/PCA/fmol_ug_glutamine_limitation_PCA.pdf", width=8, height=6)

# Gln conditions:
glc_conditions <- expr_matrix_corrected_df %>%
  dplyr::select(matches("Glc_00|Glc_06|Glc_18|Glc_29|Full_Med"))
# PCA
pca_res <- prcomp(t(glc_conditions))

pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- colnames(glc_conditions)
pca_df$condition <- sub("_[^_]+$", "", pca_df$sample)

# Percent variance
percentVar <- round(100 * summary(pca_res)$importance[2, 1:5], 1)

# Plot
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3,aes(color = condition)) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = cbPalette) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )+
  ggtitle("PCA plot of glucose limitation samples (TMT proteomics)")
#ggsave("6_results/Proteomics/PCA/fmol_ug_glucose_limitation_PCA.pdf", width=8, height=6)



