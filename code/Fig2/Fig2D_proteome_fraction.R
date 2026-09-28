## Plot proteome allocation
# Written by: Cyriel Huijer
# Date: 29 January 2026
# Last updated: 19 June 2026

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(GO.db)

# clear environment
rm(list=ls())

# Set working directory to /code/data directory
setwd(dir = "")

# Format:

# Load in abs quant data:
protein <- read.csv(file = "protein_fmol_cell.csv")
row.names(protein) <- protein$gene_name
protein <- protein %>% dplyr::select(!gene_name)
rna <- read.csv(file = "rna_amol_cell.csv")
rna <- rna %>% filter(!is.na(gene_name))
row.names(rna) <- rna$gene_name
rna <- rna %>% dplyr::select(!gene_name)

common_genes <- intersect(rownames(protein), rownames(rna))

# Formatting of conditions 
rna <- rna %>% dplyr::select(-c(Glc_00_gL_1,Glc_70_gL_1,Glc_70_gL_2,Glc_70_gL_3,Gln_00_mM_3,Gln_00_mM_2,Gln_00_mM_1)) # test
rna <- rna %>% dplyr::select(contains("Gln"),contains("Full_Med")) %>%filter(apply(., 1, var) != 0)

conditions <- rna %>% 
  colnames() %>% 
  str_remove("_\\d+$") %>% 
  unique()
protein_avg <- sapply(conditions, function(cond) {
  cols <- grep(paste0("^", cond, "_\\d+$"), colnames(protein), value = TRUE)
  rowMeans(protein[, cols, drop = FALSE])
}) %>%
  as.data.frame()

protein_long <- protein_avg %>%
  mutate(gene_name = rownames(.)) %>%
  pivot_longer(
    cols = -gene_name,
    names_to = "sample_name",
    values_to = "protein_fmol_per_cell"
  )

# Filter out low abundant proteins between cutoff
protein_long <- protein_long %>%
  filter(protein_fmol_per_cell > 1e-8)

genes <- rownames(rna)
genes2 <- rownames(protein)
genes <- union(genes, genes2)

# Allocation part:
gaf <- fread("Fig2/manual_GOSlim_top_terms.gaf", sep = "\t", header = F)
colnames(gaf) <- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID","DB_Reference","Evidence_Code","With_From","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon","Date","Assigned_By","Annotation_Extension","Gene_Product_Form_ID")

goterms <- Term(GOTERM)
go_ids <- unique(gaf$GO_ID)

go_names <- Term(GOTERM[go_ids])

# Make a table with go term names and IDs
go_table <- data.frame(
  GO_ID   = names(go_names),
  GO_name = as.character(go_names),
  stringsAsFactors = FALSE
)

gaf_with_name <- merge(gaf, go_table, by.x="GO_ID", by.y="GO_ID", all.x=TRUE)

# Only select genes that are in the dataset:
gaf_filtered <- gaf_with_name %>%
  filter(DB_Object_Symbol %in% genes)

# count occurrences per gene per GO term
counts <- gaf_filtered %>%
  group_by(DB_Object_Symbol, GO_ID) %>%
  summarise(n = n(), .groups = "drop")

# pick the most common GO term per gene
majority <- counts %>%
  group_by(DB_Object_Symbol) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup()

#ensure all genes are present, add "other" if no match
gene_assignment <- data.frame(Gene = genes) %>%
  left_join(majority, by = c("Gene" = "DB_Object_Symbol")) %>%
  mutate(GO_ID = ifelse(is.na(GO_ID), "other", GO_ID))

# Merge into gene_assignment
gene_assignment <- gene_assignment %>%
  left_join(go_table, by = "GO_ID")

# Also set genes to other
gene_assignment$GO_name[gene_assignment$GO_ID == "other"] <- "other"

## Calculate proteome fraction

# Annotate RNA with GO
protein_annotated <- protein_long %>%
  left_join(gene_assignment, by = c("gene_name" = "Gene"))
protein_annotated <- protein_annotated %>% mutate(GO_name = if_else(GO_name == "cell cycle", "other", GO_name)) # cell cycle is not shown in final figure

protein_go_summary <- protein_annotated %>%
  group_by(sample_name) %>%
  mutate(total_amol = sum(protein_fmol_per_cell, na.rm = TRUE),
         pct_of_total = (protein_fmol_per_cell / total_amol) * 100) %>%
  ungroup() %>%
  group_by(sample_name, GO_name) %>%
  summarise(total_fmol = sum(protein_fmol_per_cell, na.rm = TRUE),
            pct_of_total = sum(pct_of_total, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(sample_name = factor(sample_name,
                              levels = c("Gln_04_mM", "Gln_12_mM", "Gln_20_mM",
                                         "Gln_28_mM", "Full_Med")
  ))

protein_go_summary$GO_name <- factor(protein_go_summary$GO_name, levels=c('cell cycle', 'metabolic process', 'response to stress','translation','other'))
label_map <- c(
  "Gln_04_mM" = "Gln (0.04 mM)",
  "Gln_12_mM" = "Gln (0.12 mM)",
  "Gln_20_mM" = "Gln (0.20 mM)",
  "Gln_28_mM" = "Gln (0.28 mM)",
  "Full_Med" = "Full Media (2.5 mM)"
)
protein_go_summary
colors <- c("#D81B60","#1E88E5","#FFC107","#004D40","#920BF3")

# Fig 2D
# Plot % stacked
ggplot(protein_go_summary, aes(x = sample_name, y = pct_of_total, fill = GO_name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "% mol/mol", x = "Condition",
       fill = "GO function") +
  scale_fill_manual(values=colors)+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(labels = label_map) +  coord_flip() +
  geom_text(data = subset(protein_go_summary, pct_of_total > 2),
            aes(label = sprintf("%.1f%%", pct_of_total)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white")
# ggsave(
#   paste0("./6_results/final/fig2/", format(Sys.Date(), "%Y%m%d"), "_proteome_fraction_molmol"),
#   width = 6.5, height = 3
# )