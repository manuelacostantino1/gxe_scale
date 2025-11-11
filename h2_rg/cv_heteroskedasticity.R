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
htsk_file <- read.table("/Users/manuelacostantino/Documents/descaler_inputs/heteroskedasticity_all_phenos.txt", header = TRUE)
head(htsk_file)

htsk_file <- htsk_file[!(htsk_file$pheno=="Arm_fat-free_mass_left"),]
htsk_file <- htsk_file[!(htsk_file$pheno=="Arm_fat-free_mass_right"),]
#htsk_file <- htsk_file[!(htsk_file$pheno=="Testosterone"),]

# Merge the dfs
htsk_file$phenos <- tolower(htsk_file$pheno)
env_df$phenos <- tolower(env_df$phenotype)
merged_df <- merge(htsk_file,env_df, by= "phenos")
merged_df <- merge(merged_df,cv_df, by= "phenos")
head(merged_df)

# Make factor with low and high cv
merged_df <- merged_df %>%
  mutate(
    category = ifelse(cv < 0.2, "Low CV (<0.2)", "High CV (>0.2)"),
    category = factor(category, levels = c("Low CV (<0.2)", "High CV (>0.2)"))
  )

merged_df$hetero_log <- merged_df$log_var_female/merged_df$log_var_male
merged_df$hetero_def <- merged_df$var_female/merged_df$var_male
merged_df$yaxis <- abs(merged_df$hetero_def-merged_df$hetero_log)/merged_df$hetero_def
merged_df$phenotype <- gsub("_", " ", merged_df$phenotype)
merged_df$phenotype <- capitalize_first(merged_df$phenotype)

# Make plot
png("/Users/manuelacostantino/Desktop/cv_plot_heteroskedasticity.png",width = 10, height = 6, units = "in", res = 600)
ggplot(merged_df, aes(x = abs(coef), y = yaxis, label = phenotype)) +
  geom_point(size = 2, alpha = 0.8, color="steelblue") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black") +
  facet_wrap(~ category, ncol = 2) +
  labs(
    x = "|Main Effect of Sex|",
    y = "|% Change in Var_fem/Var_male between log and default|"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12)
  )
dev.off()


ggplot(merged_df, aes(x = cv, y = yaxis, label = phenotype)) +
  geom_point(size = 2, alpha = 0.8, color = "steelblue") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Coefficient of Variation (CV)",
    y = "|% Change in Var_fem/Var_male between log and default|"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank()
  )




