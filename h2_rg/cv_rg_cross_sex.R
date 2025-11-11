rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(tibble)
library(ggrepel)


source("/Users/manuelacostantino/Documents/descaler/local_scripts/plot_functions.R", echo = TRUE, print.eval = TRUE)
cv_df <- read.table("/Users/manuelacostantino/Documents/descaler_inputs/cvs.txt", header = TRUE)
env_df <- read.table("/Users/manuelacostantino/Documents/descaler_inputs/env_main_effects_sex.txt", header = TRUE)
rgs_file <- read.csv("/Users/manuelacostantino/Documents/descaler_inputs/rgs_cross_sex.csv", header = TRUE)
head(rgs_file)

rgs_file <- rgs_file[!(rgs_file$phenotype=="Arm_fat-free_mass_left"),]
rgs_file <- rgs_file[!(rgs_file$phenotype=="Arm_fat-free_mass_right"),]
#rgs_file <- rgs_file[!(rgs_file$pheno=="Testosterone"),]
rgs_file <- rgs_file[!(rgs_file$pheno=="WHRadjBMI_Emdin"),]
#rgs_file <- rgs_file[!(rgs_file$pheno=="C-reactive_protein"),]

# Merge the dfs
rgs_file$phenos <- tolower(rgs_file$phenotype)
env_df$phenos <- tolower(env_df$phenotype)
merged_df <- merge(rgs_file,env_df, by= "phenos")
merged_df <- merge(merged_df,cv_df, by= "phenos")
head(merged_df)

# Make factor with low and high cv
merged_df <- merged_df %>%
  mutate(
    category = ifelse(cv < 0.3, "Low CV (<0.3)", "High CV (>0.3)"),
    category = factor(category, levels = c("Low CV (<0.3)", "High CV (>0.3)"))
  )

merged_df <- na.omit(merged_df)
hist(as.numeric(merged_df$rg_se_log))
View(merged_df)
merged_df$rg_default <- as.numeric(merged_df$rg_default)
merged_df$rg_log <- as.numeric(merged_df$rg_log)
merged_df$yaxis <- abs(merged_df$rg_default-merged_df$rg_log/merged_df$rg_default)
merged_df$phenos <- gsub("_", " ", merged_df$phenos)
merged_df$phenos <- capitalize_first(merged_df$phenos)
merged_df <- merged_df %>%
  mutate(phenos = case_when(
    phenos == "Igf-1" ~ "IGF-1",
    phenos == "Hdl" ~ "HDL",
    phenos == "Rbc" ~ "RBC",
    phenos == "Whradjbmi zhu" ~ "WHRadjBMI",
    phenos == "Ea4" ~ "Education",
    phenos == "Hba1c" ~ "HbA1c",
    phenos == "Bmi" ~ "BMI",
    phenos == "Vitamin d" ~ "Vitamin D",
    TRUE ~ phenos
  ))

# Make plot
ggplot(merged_df, aes(x = abs(coef), y = yaxis, label = phenos, color=cv)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black") +
  facet_wrap(~ category, ncol = 2) +
  labs(
    x = "|Main Effect of Sex|",
    y = "|% Change in cross-sex rg between log and default|"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )



png("/Users/manuelacostantino/Desktop/cv_plot.png",width = 10, height = 6, units = "in", res = 300)
ggplot(merged_df, aes(x = as.numeric(cv), y = as.numeric(yaxis), label = phenos, color = abs(coef))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black") +
  labs(
    x = "Coefficient of Variation",
    y = "Relative change in cross-sex rg (log vs. default)",
    color = "|Main Effect of Sex|"
  ) +
  scale_color_gradient(
    low = "gray25",
    high = "#3366FF",
    guide = guide_colorbar(
      title.position = "right",   # moves title above bar
      title.theme = element_text(angle = 270),
      title.hjust = 0.5,       # centers title
      barwidth = 1,
      barheight = 24
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()

lm_model <- lm(yaxis ~ cv, data = merged_df)
p_value <- summary(lm_model)$coef["cv", "Pr(>|t|)"]
p_value


