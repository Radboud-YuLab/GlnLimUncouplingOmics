## Check protein/cell amount of protein harvest experiment. 
# Written by: Cyriel Huijer
# Last updated: 19 June 2026

library(tidyverse)
library(readxl)

# Set working directory to /code/data directory
setwd("")

# Clear environment:
rm(list=ls())

# Load in data:
data <- read_excel('cell_protein_RNA_overview.xlsx')
bars_order <- c("Gln_00_mM", "Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med", "Glc_00_gL", "Glc_06_gL", "Glc_18_gL", "Glc_29_gL", "Glc_70_gL")

# Factorize:
data <- data %>%
  mutate(condition = factor(condition, levels = bars_order))

# Keep only the relevant conditions
plot_data <- data %>%
  filter(condition %in% c("Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med"))

# Relabel conditions for prettier x-axis
plot_data$condition <- factor(plot_data$condition,
                              levels = c("Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med"),
                              labels = c("Gln (0.04 mM)", "Gln (0.12 mM)", "Gln (0.20 mM)", "Gln (0.28 mM)", "Full Media (2.5 mM)"))


plot_data_rna <- plot_data
plot_data_protein <- plot_data %>% 
  filter(!replicate == 2)

# Fig 1B
ggplot(plot_data_rna, aes(x = condition, y = ng_RNA_cell)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7, fill = "gray70") +  
  geom_point(size = 2, alpha = 0.7) +     
  labs(
    x = "Condition",
    y = "RNA/cell (ng)",
    title = ""
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(
#   paste0("../6_results/final/fig1/", format(Sys.Date(), "%Y%m%d"), "_rna_cell.pdf"),
#   width = 5, height = 5
# )

# Fig 1C
ggplot(plot_data_protein, aes(x = condition, y = ng_cell_gr)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7, fill = "gray70") +  
  geom_point(size = 2, alpha = 0.7) +     
  labs(
    x = "Condition",
    y = "Protein/cell (ng)",
    title = ""
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(
#   paste0("../6_results/final/fig1/", format(Sys.Date(), "%Y%m%d"), "_protein_cell.pdf"),
#   width = 5, height = 5
# )


































# Plot:

ggplot(data, aes(x = condition, y = ng_cell)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7) +  
  geom_point(size = 2, alpha = 0.6, aes(shape = as.factor(replicate))) +     
  labs(
    x = "Condition",
    y = "Protein per cell (ng/cell)",
    title = "Protein per cell across glutamine/glucose conditions"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("../6_results/biomass_per_condition/protein_per_cell.pdf", width=10, height=6)

# Do ANOVA to test significance:

anova_result <- aov(ug_cell ~ condition, data = data)
summary(anova_result)

Tukey <- TukeyHSD(anova_result)

Tukey_df <- as.data.frame(Tukey$condition)
Tukey_df$comparison <- rownames(Tukey$condition)


# NO significant result between replicates!


# Plot total protein harvested:

data <- data %>%
  filter(!replicate == 2)

ggplot(data, aes(x = condition, y = ug_total)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7) +  
  geom_point(size = 2, alpha = 0.6, aes(shape = as.factor(replicate))) +     
  labs(
    x = "Condition",
    y = "Total protein at harvest (ug)",
    title = "Total protein across glutamine/glucose conditions"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Also test significance:

anova_result <- aov(ug_total ~ condition, data = data)
summary(anova_result)

Tukey <- TukeyHSD(anova_result)

Tukey_df <- as.data.frame(Tukey$condition)
Tukey_df$comparison <- rownames(Tukey$condition)

plot(Tukey)

# Check protein/cell using the GR data from 20250507_CyQUANT_growth_curves.R:

ggplot(data, aes(x = condition, y = ng_cell_gr)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7) +  
  geom_point(size = 2, alpha = 0.6, aes(shape = as.factor(replicate))) +     
  labs(
    x = "condition",
    y = "protein/cell in ng",
    title = "ng protein/cell across glutamine/glucose conditions"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("../6_results/biomass_per_condition/protein_per_cell.pdf", width=10, height=6)


# Do anova:
anova_result <- aov(ng_cell_gr ~ condition, data = data)
summary(anova_result)

Tukey <- TukeyHSD(anova_result)

Tukey_df <- as.data.frame(Tukey$condition)
Tukey_df$comparison <- rownames(Tukey$condition)

plot(Tukey)


## RNA comparison

data <- read_excel('cell_protein_RNA_overview.xlsx')
bars_order <- c("Gln_00_mM", "Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med", "Glc_00_gL", "Glc_06_gL", "Glc_18_gL", "Glc_29_gL", "Glc_70_gL")
# Factorize:
data <- data %>%
  mutate(condition = factor(condition, levels = bars_order))

ggplot(data, aes(x = condition, y = ng_RNA_cell)) +
  geom_bar(stat = "summary", fun = "mean", alpha = 0.7) +  
  geom_point(size = 2, alpha = 0.6, aes(shape = as.factor(replicate))) +     
  labs(
    x = "condition",
    y = "RNA/cell in ng",
    title = "RNA/cell across glutamine/glucose conditions"
  ) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Do anova:
anova_result <- aov(ng_RNA_cell ~ condition, data = data)
summary(anova_result)

Tukey <- TukeyHSD(anova_result)

Tukey_df <- as.data.frame(Tukey$condition)
Tukey_df$comparison <- rownames(Tukey$condition)

plot(Tukey)

# I have RNA data in attomol/ug total RNA, calculate per replicate the RNA per cell (same for protein)
# Calculate through per replicate, match proteins with RNA










