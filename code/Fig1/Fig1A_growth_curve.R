## growth rate plot (Fig 1A)
# written by: Cyriel Huijer
# last updated: 19 June 2026

# Load packages
library(tidyverse)

# clear environment
rm(list=ls())

# Set working directory to /code/data directory
setwd("")

data <- read.csv("Fig1/CyQUANT_cellCounts.csv")
data <- data %>% dplyr::select(c("timepoint","Gln_04_mM", "Gln_12_mM", "Gln_20_mM", "Gln_28_mM", "Full_Med"))  

# Plot first condition with larger axis text
plot(data$timepoint, data$Gln_04_mM, type = "l", col = "#D81B60", lwd = 2,
     xlab = "Time (h)", ylab = "Cells",
     main = "",
     ylim = range(data[ ,2:6], na.rm = TRUE),
     cex.axis = 1.3,   
     cex.lab = 1.5,    
     cex.main = 1.8)   

# Add open points on top
points(data$timepoint, data$Gln_04_mM, col = "#D81B60", pch = 21, bg = "white", cex = 2, lwd = 2)

# Repeat for other conditions
lines(data$timepoint, data$Gln_12_mM, col = "#1E88E5", lwd = 2)
points(data$timepoint, data$Gln_12_mM, col = "#1E88E5", pch = 21, bg = "white", cex = 2, lwd = 2)

lines(data$timepoint, data$Gln_20_mM, col = "#FFC107", lwd = 2)
points(data$timepoint, data$Gln_20_mM, col = "#FFC107", pch = 21, bg = "white", cex = 2, lwd = 2)

lines(data$timepoint, data$Gln_28_mM, col = "#004D40", lwd = 2)
points(data$timepoint, data$Gln_28_mM, col = "#004D40", pch = 21, bg = "white", cex = 2, lwd = 2)

lines(data$timepoint, data$Full_Med, col = "#5F8247", lwd = 2)
points(data$timepoint, data$Full_Med, col = "#5F8247", pch = 21, bg = "white", cex = 2, lwd = 2)

# Legend
legend("topleft", legend = c("Gln (0.04 mM)", "Gln (0.12 mM)", "Gln (0.20 mM)", "Gln (0.28 mM)", "Full Media (2.5 mM)"),
       col = c("#D81B60","#1E88E5","#FFC107","#004D40","#5F8247"),
       lwd = 2, pch = 21, pt.bg = "white", pt.cex = 2, bty = "n", cex = 1.2)  # control legend text size
