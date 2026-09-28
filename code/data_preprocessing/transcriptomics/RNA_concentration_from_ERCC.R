## ERCC92 plots conversion
# Written by: CYriel Huijer
# Last updated: 29 June 2025

library(tidyverse)
library(purrr)
library(org.Hs.eg.db)
library(sva)
library(ggrepel)

# Set working directory to /code/data directory
setwd(dir = "")

# Load data files:
ercc_concentrations <- read.csv(file = "transcriptomics/ERCC92_concentrations.txt", sep = "\t")
counts <- read.csv(file = "transcriptomics/GRCh38.p14_ERCC92-counts.tsv", sep = "\t")
gene_lengths <- read.csv(file = "transcriptomics/GRCh38.p14_ERCC92-gene_lengths.tsv",sep = "\t")

# for translating ensembl ids to gene names:
add_genesymbol = function(genelist,decimal = TRUE){
  #Function to set ensemble names to hgnc symbol
  ##If isoform is annotated, use decimal == TRUE
  if(decimal){
    ensemble_short = as.data.frame(genelist) %>% separate(genelist,".")
  }else{
    ensemble_short = as.data.frame(genelist)
  }
  ensemble_short = AnnotationDbi::select(org.Hs.eg.db,
                                         keys=ensemble_short[,1], columns=c("SYMBOL"), keytype="ENSEMBL",multiVals = "first")
  ensemble_short = distinct(ensemble_short, ENSEMBL, .keep_all = TRUE)
  return(ensemble_short$SYMBOL)
}

counts <- merge(counts, gene_lengths,by.x = 'gene',by.y = 'gene_id')

# ERCC mix 2 was used, 6 uL of 1:1000 dilution was added to 300 ng total RNA

ercc_concentrations$attomoles_ugRNA_mix2 <- ercc_concentrations$concentration.in.Mix.2..attomoles.ul. * ((6/1000)/0.3)

# Calculation RPKM: RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 ) --> https://www.metagenomics.wiki/pdf/qc/RPKM

sample_cols <- colnames(counts)[startsWith(colnames(counts), "X")]
rpkm <- counts


## Instead of running this, open the RPKM file (it takes a few mins to run otherwise):
rpkm <- read.csv("transcriptomics/outputs/RPKM_values.csv")
# for (sample in sample_cols) {
#   print(sample)
#   for (i in 1:nrow(counts)){
#     rpkm[[sample]][i] <- counts[[sample]][i] / ( counts$length[i]/1000 * sum(counts[[sample]])/1e6) 
#   } 
# }

#write.csv(rpkm, "transcriptomics/outputs/RPKM_values.csv",row.names = F)

rpkm_ercc <- rpkm %>%
  filter(grepl("ERCC", gene))

# Identify the expression columns (excluding 'gene' and 'length')
expr_cols <- setdiff(names(rpkm_ercc), c("gene", "length"))

# Apply filter: set values < 1 or > 10000 to NA
rpkm_ercc[expr_cols] <- lapply(rpkm_ercc[expr_cols], function(col) {
  col[col < 1] <- NA
  return(col)
})

rpkm_ercc_conc <- merge(rpkm_ercc, ercc_concentrations,by.x = 'gene',by.y = 'ERCC.ID')

# Generate plots
plots <- map(sample_cols, function(sample_col) {
  
  df <- rpkm_ercc_conc %>%
    filter(.data[[sample_col]] > 0) %>%
    mutate(log10_x = log10(attomoles_ugRNA_mix2),
           log10_y = log10(.data[[sample_col]]))
  
  # Calculate R-squared
  r_squared <- round(summary(lm(log10_y ~ log10_x, data = df))$r.squared, 3)
  n <- nrow(df)
  
  # Generate plot
  ggplot(df, aes(x = log10_x, y = log10_y)) +
    geom_point(color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_text(aes(x = min(log10_x, na.rm = TRUE), 
                  y = max(log10_y, na.rm = TRUE), 
                  label = paste0("R² = ", r_squared,"\nn = ",n)),
              hjust = 0, vjust = 1, size = 5) +
    labs(title = paste0("ERCC92 spike-in (RPMK): ", sample_col),
         x = "log10(attomoles / ug RNA)",
         y = paste0("log10(RPKM): ", sample_col)) +
    theme_bw(base_size = 15)
})

# give names to plots
names(plots) <- sample_cols

# display all plots
for (p in plots) {
  print(p)
}

# plot all in one figure:

# reshape to long format for plot
rpkm_long <- rpkm_ercc_conc %>%
  dplyr::select(gene, attomoles_ugRNA_mix2, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "rpkm") %>%
  filter(rpkm > 0, attomoles_ugRNA_mix2 > 0) %>%
  mutate(log10_rpkm = log10(rpkm),
         log10_attomoles = log10(attomoles_ugRNA_mix2))

# plot with all samples overlaid
ggplot(rpkm_long, aes(x = log10_attomoles, y = log10_rpkm, color = sample)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "RPKM vs Spike-in Abundance Across Samples",
       x = "log10(attomoles / ug RNA)",
       y = "log10(RPKM)",
       color = "Sample") +
  theme_bw(base_size = 15)

## Calculation of transcript levels in each condition: #########################
rpkm_without_ercc <- rpkm %>%
  filter(!grepl("ERCC", gene))
attomol_ug <- rpkm_without_ercc

# Identify the expression columns (excluding 'gene' and 'length')
expr_cols <- setdiff(names(attomol_ug), c("gene", "length"))

# Apply filter: set values < 1 or > 10000 to NA
attomol_ug[expr_cols] <- lapply(attomol_ug[expr_cols], function(col) {
  col[col < 1 | col > 10000] <- NA
  return(col)
})


# Get attomol/ug RNA using lm for the ERCC transcripts

for (sample in sample_cols){
  print(sample)
  curve_calc <- rpkm_ercc_conc[,c('gene',sample,'attomoles_ugRNA_mix2')]
  colnames(curve_calc) <- c('gene','rpkm','attomoles_ugRNA_mix2')
  curve_calc <- curve_calc %>% filter(rpkm > 0)
  a <- coef(lm(log10(curve_calc$rpkm) ~ log10(curve_calc$attomoles_ugRNA_mix2)))[2]
  b <- coef(lm(log10(curve_calc$rpkm) ~ log10(curve_calc$attomoles_ugRNA_mix2)))[1]
  attomol_ug[[paste0(sample, "_attomoles_ugRNA")]] <- 10^((log10(attomol_ug[[sample]]) - b) / a)
}

attomol_ug <- attomol_ug %>%
  dplyr::select(gene, ends_with("_attomoles_ugRNA"))
attomol_ug <- attomol_ug[!apply(attomol_ug[, paste0(sample_cols,"_attomoles_ugRNA")], 1, function(x) all(is.na(x))), ]



attomol_ug$gene_symbol <- add_genesymbol(attomol_ug$gene)


# Exclude non-numeric columns: keep only the sample columns

sample_data <- attomol_ug[, 2:34]  # columns 2 to 34 are numeric sample columns
row.names(sample_data) <- attomol_ug$gene
sample_data <- sample_data[complete.cases(sample_data), ]

sample_data_t <- t(sample_data)

# Perform PCA
pca_res <- prcomp(sample_data_t, center = TRUE, scale. = TRUE)

# Create a data frame for plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$sample <- rownames(pca_df)  # sample names from transposed data

ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) +
  labs(title = "PCA of Samples (Absolute RNA Quantities)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")) +
  theme_bw(base_size = 15)

## Do batch correction using ComBat 

colnames(attomol_ug) <- sub("^[^_]+_", "", colnames(attomol_ug))
attomol_ug <- attomol_ug[complete.cases(attomol_ug), ]

#Create batch vector (one value per column/sample)
batch_vector <- as.character(c(
  rep(1, 11),  
  rep(2, 11),
  rep(3, 11)
))

base_conditions <- c(
  "Gln_00_mM", "Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med",
  "Glc_00_gL", "Glc_06_gL", "Glc_18_gL", "Glc_29_gL", "Glc_70_gL"
)
conditions <- rep(base_conditions, times = 3)
condition_batch = paste0(conditions,"_",batch_vector)
batch <- rep(1:3, each = 11)


# drop unnecessary columns 
rownames(attomol_ug) <- attomol_ug$gene
colnames(attomol_ug)[2:34] <- condition_batch
attomol_ug <- attomol_ug %>% dplyr::select(-gene,-symbol)
attomol_ug <- attomol_ug[, !grepl("Glc_06_gL_1", colnames(attomol_ug))] # From QC, this sample is a clear outlier, also see previous PCA, so remove it before batch correction
#attomol_ug <- attomol_ug[, !colnames(attomol_ug) %in% c("Gln_00_mM_1","Gln_04_mM_1","Glc_06_gL_1","Glc_18_gL_1","Glc_18_gL_2","Glc_18_gL_3","Glc_70_gL_1","Glc_70_gL_2","Glc_70_gL_3","Gln_12_mM_1","Gln_12_mM_2","Gln_12_mM_3","Gln_28_mM_1","Gln_28_mM_2","Gln_28_mM_3")]
attomol_ug_matrix <- as.matrix(attomol_ug)

log2_expr <- log2(attomol_ug_matrix + 1) # Log transform data: https://www.biostars.org/p/219320/

# Set condition vector as mod
mod <- model.matrix(~ as.factor((sub("_[^_]+$", "", colnames(log2_expr)))))
expr_matrix_corrected <- ComBat(dat = log2_expr, batch = sub(".*_", "", colnames(log2_expr)), mod = mod)

# Perform PCA
expr_matrix_filtered <- expr_matrix_corrected[apply(expr_matrix_corrected, 1, function(x) var(x) != 0), ]
pca_result <- prcomp(t(expr_matrix_filtered), scale. = TRUE)

# Prepare PCA results for plotting
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Batch = sub(".*_", "", colnames(attomol_ug_matrix)),
  Condition = factor(sub("_[^_]+$", "", colnames(attomol_ug_matrix))),  # Add sample names
  Sample = paste0(factor(sub("_[^_]+$", "", colnames(attomol_ug_matrix))),"_",sub(".*_", "", colnames(attomol_ug_matrix)))
)

# Plot using ggplot2 with labels and explained variance on axes
explained_variance <- summary(pca_result)$importance[2, ] * 100  # Percent of variance explained

# Use colorblind-friendly palette
cbPalette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999", "#000000", "#AAAA00", "#44AA99"
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition,fill = Condition)) +
  geom_point(size = 3,aes(fill = Condition)) +
  theme_bw(base_size = 15) +
  labs(
    title = "PCA After ComBat Correction",
    x = paste("PC1 -", round(explained_variance[1], 2), "% Explained Variance"),
    y = paste("PC2 -", round(explained_variance[2], 2), "% Explained Variance")
  ) +
  scale_color_manual(values = cbPalette) +
  theme(legend.position = "right") +
  scale_shape_manual(values=1:nlevels(pca_df$Condition))+
  geom_text_repel(aes(label = Sample), size = 4, max.overlaps = Inf)

## Select only Gln:

# Gln conditions:
expr_matrix_filtered <- as.data.frame(expr_matrix_filtered)
gln_conditions <- expr_matrix_filtered %>%
  dplyr::select(matches("Gln|Full_Med"))
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
  ggtitle("PCA plot of glutamine limitation samples (RNAseq)")
#ggsave("6_results/RNAseq/PCA/amol_ug_glutamine_limitation_PCA.pdf", width=8, height=6)

# Save batch corrected data
expr_matrix_filtered_amol_ug <- 2^expr_matrix_filtered - 1
expr_matrix_filtered_amol_ug$gene_name <- rownames(expr_matrix_filtered_amol_ug)
expr_matrix_filtered_amol_ug$gene_name <- add_genesymbol(expr_matrix_filtered_amol_ug$gene_name)
expr_matrix_filtered_amol_ug <- expr_matrix_filtered_amol_ug[rowSums(expr_matrix_filtered_amol_ug < 0) == 0, ] # remove negative values after batch correction
#write.csv(expr_matrix_filtered_amol_ug, "transcriptomics/outputs/amol_per_ug_RNA.csv",row.names = F)
