# UPS2 analysis 
# written by: Cyriel Huijer
# Last updated: 06 May 2025
# Data from: CH00032 (UPS2 run)
# Filename: 
# Control run (no UPS2) = 20240802_Exploris1_SA_CH_MCF10A_CH32_GlnGlc_control_15ug_50percent.raw
# UPS2 run = 20240802_Exploris1_SA_CH_MCF10A_CH32_GlnGlc_UPS2_18ug_50percent.raw

library(tidyverse)
library(ggpmisc)

# Set working directory to /code/data directory
setwd("")

# This file is the control file (not spiked with UPS2)
#raw_data <- read.csv("proteomics/MQ_no_UPS2.txt",sep = "\t")

# This file is the MCF10A proteome with UPS2: 
raw_data <- read.csv(file = "proteomics/MQ_UPS2.txt",sep = "\t")
ups2_meta <- read.csv("proteomics/UPS_proteome_set.csv",header = TRUE)

filtered_data <- raw_data[raw_data$Only.identified.by.site != "+" &
                            raw_data$Reverse != "+" &
                            raw_data$Potential.contaminant != "+", ]
filtered_data_with_ups <- filtered_data %>%
  filter(grepl("ups", Protein.IDs)&
           !grepl("^G", Protein.IDs))%>%
  mutate(Protein.IDs_6 = substr(Protein.IDs, 1, 6))  # Extract the first 6 characters

merged_data <- merge(filtered_data_with_ups, ups2_meta, 
                     by.x = "Protein.IDs_6", 
                     by.y = "UniProt.Accession.Number", 
                     all.x = TRUE)  # Keep all rows from filtered_data_with_ups_no_G

merged_data <- merged_data %>%
  mutate(
    iBAQ_log10 = log10(iBAQ + 1), # Log10 transformation of iBAQ
    UPS2_log10 = log10(UPS2.Amount..fmol.*(3/10.6))  # Log10 transformation of UPS2Amount (correct for protein)
  )

filtered_merged_data <- merged_data %>%
  filter(Protein.IDs_6 != "P08758" &
           Protein.IDs_6 != "P02787"&
           UPS2_log10 > 1 )

x = filtered_merged_data$UPS2_log10
y = filtered_merged_data$iBAQ_log10
formula <- y ~ x
ggplot(filtered_merged_data, aes(x = UPS2_log10, y = iBAQ_log10)) +
  geom_point(alpha = 0.6,size = 3) +  
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add a linear regression line
  labs(
    title = "UPS2 Spike-in iBAQ signal",
    x = "UPS2 (log10-transformed)",
    y = "iBAQ (log10-transformed)"
  ) +
  stat_poly_eq(use_label(c("eq","adj.R2","P")),formula=formula)+
  theme_bw(base_size = 15) 

model <- lm(iBAQ_log10~UPS2_log10, data=filtered_merged_data)
summary(model)# Adjusted R2 = 0.72
model$coefficients["(Intercept)"]


## CONTINUE HERE!!!

# Formula for calculating the protein amounts for other proteins: 
# iBAQ = 6.74 + 0.765[femtomoles] 
# To retrieve the femtomoles per protein, calculate as follows:
# [femtomoles] = (iBAQ - 6.74)/0.765

filtered_data_without_ups <- filtered_data %>%
  filter(!grepl("ups", Protein.IDs)) %>%
  mutate(Protein.IDs_6 = substr(Protein.IDs, 1, 6))

## CALCULATE ABSOLUTE VALUES FROM THIS DATA

filtered_data_without_ups <- filtered_data_without_ups %>%
  mutate(
    iBAQ_log10 = log10(iBAQ + 1),  # Log10 transformation of iBAQ
  )
filtered_data_without_ups$femtomoles <- 10^((filtered_data_without_ups$iBAQ_log10 - model$coefficients["(Intercept)"])/model$coefficients["UPS2_log10"])
filtered_data_without_ups <- filtered_data_without_ups[!grepl("ups", filtered_data$Protein.IDs), ]

filtered_data_without_ups$Accession <- gsub(".*\\|(.*)\\|.*", "\\1", filtered_data_without_ups$Fasta.headers)


quant_results <- filtered_data_without_ups[, c("Gene.names", "femtomoles")]
write.csv(quant_results, file = "proteomics/output/UPS2_quantification_results.csv", row.names = FALSE)

