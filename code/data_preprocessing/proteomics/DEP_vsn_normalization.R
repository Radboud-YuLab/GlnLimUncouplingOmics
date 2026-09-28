## DEP vsn normalization, outputs normalized quantities
# Written by: Cyriel Huijer
# Date: 29 April 2025

# It was found there is a batch effect between different TMT runs; therefore, 
# after correction of the pooled sample, the ComBat algorithm was applied for batch
# correction. 

# Set working directory to /code/data directory
setwd("C:/Users/chuijer/Documents/nextcloud_backup/PhD/Projects/Project3/manuscript/code/data")

# Load in libraries:
library(tidyverse) # tidyverse_2.0.0
library(DEP) # DEP_1.20.0

## ---------------------------------------------------------------------------## 
# All samples together:
## ---------------------------------------------------------------------------##

# Replicate 1:

# Load in data:

data_1 <- read_excel("proteomics/PD_output/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep1_tmt_50ug_50percent.xlsx")
data_1 <- data_1 %>% filter(Contaminant == FALSE)

data_1$protein_name <- sub(" GN=.*", "", data_1$Description)
data_1 <- make_unique(data_1, "Accession", "protein_name")

# Add sample information and create Summarized Experiment (SE) object

sample_info <- data.frame(
  label = c("TMT.pooled","TMT.Gln_04_mM", "TMT.Gln_12_mM", "TMT.Gln_20_mM", "TMT.Gln_28_mM", "TMT.Full_Med","TMT.Glc_00_gL","TMT.Glc_06_gL", "TMT.Glc_18_gL", "TMT.Glc_29_gL", "TMT.Glc_70_gL"),  # Add all your sample names
  condition = c("pooled","Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med","Glc_00_gL","Glc_06_gL", "Glc_18_gL", "Glc_29_gL", "Glc_70_gL"),              # Define conditions
  replicate = c(1,1,1,1,1,1,1,1,1,1,1) 
)

# Select columns

colnames(data_1)[36:46] <- sample_info$label

se <- make_se(data_1, columns = c(36:46), expdesign = sample_info)
TMT_columns <- grep("TMT.", colnames(data_1)) # Grep TMT columns
data_se_parsed <- make_se_parse(data_1, TMT_columns)

data_filt <- filter_missval(se, thr = 0)
data_norm <- normalize_vsn(data_filt)
data_norm_1 <- data_norm

plot_normalization(se, data_norm)
plot_pca(data_norm, x = 1, y = 2, n = 500, point_size = 4,label_size = 4,label = TRUE)

data_corr <- data.frame(data_norm@assays@data@listData[[1]]) # Gives dataframe with normalized values
#write.table(data_corr,"proteomics/output/rep1_vsn_normalized_correctlabels.csv",row.names = TRUE, sep = ',')

rm(list=ls())

# Replicate 2:

# Laod in data:

data_2 <- read_excel("proteomics/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep2_tmt_50ug_50percent.xlsx")
data_2 <- data_2 %>% filter(Contaminant == FALSE)

data_2$protein_name <- sub(" GN=.*", "", data_2$Description)
data_2 <- make_unique(data_2, "Accession", "protein_name")

sample_info <- data.frame(
  label = c("TMT.pooled","TMT.Gln_00_mM", "TMT.Gln_04_mM", "TMT.Gln_12_mM", "TMT.Gln_20_mM", "TMT.Gln_28_mM", "TMT.Full_Med","TMT.Glc_00_gL","TMT.Glc_06_gL", "TMT.Glc_18_gL", "TMT.Glc_29_gL", "TMT.Glc_70_gL"),  # Add all your sample names
  condition = c("pooled", "Gln_00_mM", "Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med","Glc_00_gL", "Glc_06_gL", "Glc_18_gL", "Glc_29_gL", "Glc_70_gL"),              # Define conditions
  replicate = c(2,2,2,2,2,2,2,2,2,2,2,2) 
)

colnames(data_2)[36:47] <- sample_info$label


se <- make_se(data_2, columns = c(36:47), expdesign = sample_info)
TMT_columns <- grep("TMT.", colnames(data_2)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_2, TMT_columns)

data_filt <- filter_missval(se, thr = 0)

data_norm <- normalize_vsn(data_filt)
data_norm_2 <- data_norm

plot_normalization(se, data_norm)
plot_pca(data_norm, x = 1, y = 2, n = 500, point_size = 4,label_size = 4,label =TRUE)

data_corr <- data.frame(data_norm@assays@data@listData[[1]]) # Gives dataframe with normalized values
#write.table(data_corr,"proteomics/output/vsn_normalized_data/rep2_vsn_normalized_correctlabels.csv",row.names = TRUE, sep = ',')

rm(list=ls())


# Replicate 3:

data_3 <- read_excel("proteomics/20240801_Exploris1_SA_CH_MCF10A_GlnGlc_rep3_tmt_50ug_50percent.xlsx")
data_3 <- data_3 %>% filter(Contaminant == FALSE)

data_3$protein_name <- sub(" GN=.*", "", data_3$Description)
data_3 <- make_unique(data_3, "Accession", "protein_name")

sample_info <- data.frame(
  label = c("TMT.Glc_70_gL","TMT.pooled","TMT.Gln_00_mM", "TMT.Gln_04_mM", "TMT.Gln_12_mM", "TMT.Gln_20_mM", "TMT.Gln_28_mM", "TMT.Full_Med","TMT.Glc_00_gL","TMT.Glc_06_gL", "TMT.Glc_18_gL", "TMT.Glc_29_gL"),  # Add all your sample names
  condition = c("Glc_70_gL","pooled", "Gln_00_mM", "Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med","Glc_00_gL", "Glc_06_gL", "Glc_18_gL", "Glc_29_gL"),              # Define conditions
  replicate = c(3,3,3,3,3,3,3,3,3,3,3,3) 
)

colnames(data_3)[c(36:46,48)] <- sample_info$label

se <- make_se(data_3, columns = as.integer(c(36:46,48)), expdesign = sample_info)
TMT_columns <- grep("TMT.", colnames(data_3)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_3, TMT_columns)

data_filt <- filter_missval(se, thr = 0)
data_norm <- normalize_vsn(data_filt)

plot_normalization(se, data_norm)

plot_pca(data_norm, x = 1, y = 2, n = 500, point_size = 4,label_size = 4,label =TRUE)

data_corr <- data.frame(data_norm@assays@data@listData[[1]]) # Gives dataframe with normalized values
#write.table(data_corr,"proteomics/output/rep3_vsn_normalized_correctlabels.csv",row.names = TRUE, sep = ',')








