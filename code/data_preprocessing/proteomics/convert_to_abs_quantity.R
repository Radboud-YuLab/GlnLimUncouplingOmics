## Calculates absolute protein quantity:
# Written by: Cyriel Huijer
# Last updated 26 May 2025

library(tidyverse)
library(gridExtra) 
library(readxl)
library(ggbreak)


rm(list=ls())

# Set working directory to /code/data directory
setwd("")

# Load in data:

#data <- read.csv("5_processedData/Proteomics/vsn_ComBat_correction/vsn_combat_corrected_log2_pool_ratios_allsamples_correctlabels.csv")
data <- read.csv("proteomics/output/vsn_combat_corrected_log2_pool_ratios_correctlabels.csv")

abs_quant <- read.csv("proteomics/output/UPS2_quantification_results.csv")
colnames(abs_quant) <- c("gene_name","femtomoles")

# Load original dataframes for GN identifiers:
data_1 <- read_excel("proteomics/PD_output/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep1_tmt_50ug_50percent.xlsx")
data_2 <- read_excel("proteomics/PD_output/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep2_tmt_50ug_50percent.xlsx")
data_3 <- read_excel("proteomics/PD_output/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep3_tmt_50ug_50percent.xlsx")

combined_data <- rbind(
  data_1[, c("Accession", "Description","MW [kDa]")],
  data_2[, c("Accession", "Description","MW [kDa]")],
  data_3[, c("Accession", "Description","MW [kDa]")]
)

unique_data <- unique(combined_data)

# Merge DFs 

data_protein_GN <- merge(data, unique_data, by = 'Accession')
data_protein_GN$gene_name <- sub(".*GN=([^ ]+) PE.*", "\\1", data_protein_GN$Description)

# Merge with absolute quant data:

merged_abs_quant <- merge(data_protein_GN,abs_quant, by = "gene_name")

# calculate absolute abundance:
condition_cols <- grep("^(Gln|Glc|Full)_", names(merged_abs_quant), value = TRUE)

# Create new columns for actual femtomoles per condition, normalized by 15 ug (total sample load)
for (col in condition_cols) {
  new_col <- paste0(col, "_fmol")
  merged_abs_quant[[new_col]] <- (merged_abs_quant$femtomoles * 2^merged_abs_quant[[col]]) / 15
}

# format absolute quantities table:
protein_fmol_ug <- merged_abs_quant %>% dplyr::select(gene_name,dplyr::ends_with("_fmol"))
#write.csv(protein_fmol_ug,"proteomics/output/fmol_ug_protein.csv",row.names = F)
