rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)

# Read in arrays
source("/Users/manuelacostantino/Documents/descaler/local_scripts/plot_functions.R", echo = TRUE, print.eval = TRUE)
load("/Users/manuelacostantino/Documents/descaler_inputs/heatmap_array.Rdata")

df <- data.frame(environment = character(), independent_count = numeric(), dependent_count = numeric(), null_count=numeric())

short_arrays <- make_3_arrays_v2(heatmap_array, threshold, lambdashort)
killable_array <- short_arrays$killable
null_array <- short_arrays$null
unkillable_array <- short_arrays$unkillable

for (env in c("sex", "age", "Smoking_status", "Alcohol_intake_frequency", "Statins", "EA4")){
  
  killable01 <- collapse_interactions(killable_array, env)
  unkillable01 <- collapse_interactions(unkillable_array, env)
  null01 <- collapse_interactions(null_array, env)

  killable_count   <- sum(killable01 == 1, na.rm = TRUE)
  unkillable_count <- sum(unkillable01 == 1, na.rm = TRUE)
  null_count       <- sum(null01 == 1, na.rm = TRUE)
  print(env)
  print(killable_count)
  print(unkillable_count)
  print(null_count)

  df <- rbind(df, c(env, unkillable_count, killable_count, null_count))
}

colnames(df) <- c("Environment", "Scale-independent\n (default+log)", "Scale-dependent\n (default only)", "null_count")
df_long <- melt(df,
                id.vars = "Environment",
                measure.vars = c("Scale-independent\n (default+log)", 
                                 "Scale-dependent\n (default only)"),
                variable.name = "Type",
                value.name = "Count")

df_long$Count <- as.numeric(df_long$Count)

# order by scale dependent count
scale_dep_counts <- df_long %>% 
  filter(Type == "Scale-dependent\n (default only)") %>%
  select(Environment, Count)

df_long$Environment <- factor(df_long$Environment,
                              levels = scale_dep_counts$Environment[order(scale_dep_counts$Count, decreasing = TRUE)])

levels(df_long$Environment)[levels(df_long$Environment) == "EA4"] <- "Educational Attainment"
levels(df_long$Environment)[levels(df_long$Environment) == "Statins"] <- "Statin Use"
levels(df_long$Environment)[levels(df_long$Environment) == "sex"] <- "Sex"
levels(df_long$Environment)[levels(df_long$Environment) == "age"] <- "Age"
levels(df_long$Environment)[levels(df_long$Environment) == "Smoking_status"] <- "Smoking Status"
levels(df_long$Environment)[levels(df_long$Environment) == "Alcohol_intake_frequency"] <- "Alcohol Intake Frequency"

# Barplot
png(paste0("/Users/manuelacostantino/Desktop/cross_env.png"), width = 6, height = 8, units = "in", res = 600)
ggplot(df_long, aes(x = Environment, y = Count, fill = Type)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Environment",
    y = "Counts",
    fill = "Type"
  ) +
  scale_fill_manual(values = c("Scale-independent\n (default+log)" = "#660099",
                               "Scale-dependent\n (default only)" = "#9999FF")) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
    axis.text.y = element_text(color="black"),
    #axis.title.x = element_text(hjust = 0.1, margin = margin(t = -10)),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white")
  )
dev.off()
