rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

# Read in arrays
source("/Users/manuelacostantino/Documents/descaler/local_scripts/plot_functions.R", echo = TRUE, print.eval = TRUE)
load("/Users/manuelacostantino/Documents/descaler_inputs/heatmap_array.Rdata")

short_arrays <- make_3_arrays_v2(heatmap_array, threshold, lambdashort)
killable_array <- short_arrays$killable
null_array <- short_arrays$null
unkillable_array <- short_arrays$unkillable

envs <- c("sex")

killable01 <- collapse_interactions(killable_array, envs)
unkillable01 <- collapse_interactions(unkillable_array, envs)
unkillable_long <- make_long(unkillable01)
killable_long <- make_long(killable01)

unkillable_long$Category <- factor(unkillable_long$Category,levels = c(0, 1), labels = c("Null", "Scale-independent\n (default + log)"))
killable_long$Category <- factor(killable_long$Category,levels = c(1, 0), labels = c( "Scale-dependent\n (default only)", " "))

killable_long <- killable_long[!(killable_long$phenotype=="WHRadjBMI Emdin"),]
unkillable_long <- unkillable_long[!(unkillable_long$phenotype=="WHRadjBMI Emdin"),]

# Switch alphabetical order
levels_all <- sort(unique(c(as.character(unkillable_long$phenotype),
                            as.character(killable_long$phenotype))))
rev_levels <- rev(levels_all)
unkillable_long$phenotype <- factor(as.character(unkillable_long$phenotype),
                                    levels = rev_levels)
killable_long$phenotype   <- factor(as.character(killable_long$phenotype),
                                    levels = rev_levels)
# Make vertical plot
png(paste0("/Users/manuelacostantino/Desktop/heatmap_", envs[1], ".png"), width = 10, height = 8, units = "in", res = 600)
ggplot() + 
  geom_tile(data = unkillable_long,
            aes(x = phenotype, y = interaction, fill = Category, color = Category)) +
  geom_tile(data = killable_long,
            aes(x = phenotype, y = interaction, fill = Category, color = Category)) +
  scale_fill_manual(
    values = c("Null" = "grey95",
               "Scale-dependent\n (default only)" = "#9999FF",
               "Scale-independent\n (default + log)" = "#660099",
               " " = "transparent"),
    guide = guide_legend(
      nrow = 1,
      keywidth = unit(2, "cm"),
      keyheight = unit(0.7, "cm"),
      byrow = TRUE,
      override.aes = list(color = NA)  
    )
  ) +
  scale_color_manual(
    values = c("Null" = "black",
               "Scale-dependent\n (default only)" = "black",
               "Scale-independent\n (default + log)" = "black",
               " " = "transparent"),
    guide = "none"
  ) +
  labs(x = "Polygenic Scores", y = "Phenotypes",
       title = paste("Gene x", capitalize_first(envs[1]), "Interactions")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.x = unit(0.8, "cm"),
    legend.box.margin = margin(4, 4, 4, 4),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
dev.off()
