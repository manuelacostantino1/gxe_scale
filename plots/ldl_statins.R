library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)


LDL <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/LDL/LDL674178.pheno", header=TRUE)
statins <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Statins/Statins674178.pheno", header=TRUE)
pgs <- read.table("/gpfs/data/ukb-share/dahl/manuela/ukb_prs_downloading/polygenic_scores/sPRS_LDL/sPRS_LDL678224.pheno", header=TRUE)

df <- merge(LDL, statins, by=c("FID","IID"))
df <- merge(df, pgs, by=c("FID","IID"))

# Rename columns for clarity
colnames(df) <- c("FID", "IID", "LDL", "statins", "pgs")

# Convert statins to factor and drop NAs
df$statins <- as.factor(df$statins)
df <- na.omit(df)

head(df)
head(df)


png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/ldl_statins.png")
ggplot(df, aes(x = pgs, y = LDL, color = statins)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +   # regression lines
  scale_color_manual(values = c("0" = "#FFCCFF", "1" = "#6600CC"),
                     labels = c("No Statins", "On Statins")) +
  labs(
    x = "Polygenic Score (LDL)",
    y = "Measured LDL (mmol/L)",
    color = "Statins",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),   # remove grid
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

dev.off()