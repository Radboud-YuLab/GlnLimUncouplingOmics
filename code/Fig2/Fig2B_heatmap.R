## Plot heatmap of different conditions (Gln lim)
# Written by: Cyriel Huijer
# Date: 9 June 2025
# Last updated: 19 June 2026

library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(clusterProfiler)
library(org.Hs.eg.db)

# Set working directory to /code/data directory
setwd("")

# clear environment
rm(list=ls())

# Load dataframes:

# Load df that includes all genes available for transcriptomics (without filter for matching protein/RNA)
subset_gln <- read.csv('Fig2/all_RNA_heatmap.csv')
subset_gln <- subset_gln %>% filter(sample_name %in% c("Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med"))
top_genes <- read.csv('Fig2/pc2_genes_heatmap.csv')
# Filter only gene names that are in PC2
subset_gln <- subset_gln %>%
  filter(gene_name %in% top_genes$gene)


# Pivot to wide format: rows = genes, columns = samples
heatmap_data <- subset_gln %>%
  pivot_wider(names_from = sample_name, values_from = rna_amol_per_cell) %>%
  column_to_rownames("gene_name")

# Log-transform
heatmap_data_log <- log10(heatmap_data + 1e-6)

# Scale rows
heatmap_scaled <- t(scale(t(heatmap_data_log)))
heatmap_scaled_df <- as.data.frame(heatmap_scaled)
heatmap_scaled_df$gene <- row.names(heatmap_scaled)

# sample order in heatmap
desired_order <- c("Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med")

# Fig 2B
# Plot heatmap
Heatmap(
  heatmap_scaled[, desired_order, drop = FALSE],
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = desired_order,
  show_row_names = TRUE,
  column_title = "Glutamine limitation",
  row_title = "Top PC2 genes",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

## Determination of nr of clusters + GO term analysis

# silhouette plot to detemrine max number of clusters
set.seed(123)
fviz_nbclust(heatmap_scaled, kmeans, method = "silhouette") + 
  labs(title = "Silhouette method for optimal k")

set.seed(123)
best_k <- 2
km_res <- kmeans(heatmap_scaled, centers = best_k, nstart = 25)

# Assign clusters to each gene
gene_clusters <- data.frame(
  gene = rownames(heatmap_scaled),
  cluster = km_res$cluster
)

head(gene_clusters)

fviz_cluster(km_res, data = heatmap_scaled,
             geom = "point", ellipse.type = "convex",
             show.clust.cent = TRUE) +
  labs(title = paste("K-means clustering (k =", best_k, ")"))

# Cut tree into 2 clusters
row_clusters <- cutree(hclust(dist(heatmap_scaled)), k = 2)

# Assign cluster labels back to genes
cluster_assignments <- data.frame(
  gene = rownames(heatmap_scaled),
  cluster = row_clusters
)


row_ha <- rowAnnotation(
  cluster = as.factor(row_clusters),
  col = list(cluster = c("1" = "red", "2" = "blue"))
)

Heatmap(
  heatmap_scaled[, desired_order, drop = FALSE],
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_order = desired_order,
  show_row_names = TRUE,
  row_title = "Top PC2 genes",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
) + row_ha



genes_by_cluster <- split(gene_clusters$gene, gene_clusters$cluster)


go_results <- list()

for (cl in names(genes_by_cluster)) {
  # Convert gene symbols to Entrez IDs
  gene_symbols <- genes_by_cluster[[cl]]
  gene_ids <- bitr(gene_symbols, fromType="SYMBOL",
                   toType="ENTREZID",
                   OrgDb="org.Hs.eg.db")
  
  # GO enrichment
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont= "BP",# Biological Process
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05)
  
  go_results[[cl]] <- ego
}

plot_cluster_go <- function(go_res, cluster_id, top_n = 20) {
  df <- as.data.frame(go_res)
  if (nrow(df) == 0) return(NULL)  # skip if no results
  
  # Prepare top GO terms
  top_df <- df %>%
    slice_head(n = top_n) %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Description = factor(Description, 
                           levels = Description[order(GeneRatio_num, decreasing = FALSE)])
    )
  
  # Bubble plot
  p <- ggplot(top_df, aes(x = Description, y = GeneRatio_num)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
    scale_size_continuous(name = "Gene Count") +
    coord_flip() +
    labs(
      x = "",
      y = "Gene Ratio",
      title = paste("Cluster", cluster_id, "Genes"),
      subtitle = "GO Biological Process (BP)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text()
    )
  
  # Return both
  return(list(
    df = df, # full results
    top_df = top_df, # top n results
    plot = p # bubble plot
  ))
}


# plot cluster 1
res_cluster1 <- plot_cluster_go(go_results[["1"]], cluster_id = 1, top_n = 20)
res_cluster1$plot

#write.csv(x = res_cluster1$df,file = "../../supplementary_tables/TableS1.csv",row.names = F)

#plot cluster 2
res_cluster2 <- plot_cluster_go(go_results[["2"]], cluster_id = 2, top_n = 20)
res_cluster2$plot

#write.csv(x = res_cluster2$df,file = "../../supplementary_tables/TableS2.csv",row.names = F)