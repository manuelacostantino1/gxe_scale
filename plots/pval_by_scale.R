rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggrepel)

FULL_RANGE <- FALSE
pheno <- "Height"
#envs <- c("sex","age", "Alcohol_intake_frequency", "Statins", "EA4", "Smoking_status")
envs <-c("sex")

# Read in arrays
source("/Users/manuelacostantino/Desktop/descaler_analysis/plot_functions.R")
load("/Users/manuelacostantino/Documents/descaler_inputs/heatmap_array.Rdata")
threshold <- 0.05/(length(prs_list))
# Make 2 versions of the arrays
pval_array <- heatmap_array
long_arrays <- make_3_arrays(pval_array, threshold, lambdaseq)
short_arrays <- make_3_arrays(pval_array, threshold, lambdashort)


# Create empty dataframe with lambda columns and a context column
if (FULL_RANGE == FALSE) {
  lambda_indices <- match(dimnames(pval_array)[[4]], as.character(lambdashort))
  lambda_indices <- which(!is.na(lambda_indices))
  pval_array <- pval_array[,,,lambda_indices]
  dimnames(pval_array)[[4]] <- as.character(lambdashort)
  lambdas <- lambdashort
} else {
  lambdas <- lambdaseq
}
df <- data.frame(matrix(nrow = 0, ncol = length(lambdas) + 2))
colnames(df) <- c(as.character(lambdas), "context", "row_id")

row_id_counter <- 1
for (env in envs) {
  subdf <- as.data.frame(pval_array[, pheno, env, ])
  colnames(subdf) <- as.character(lambdas)
  subdf$context <- env
  subdf$row_id <- seq(row_id_counter, length.out = nrow(subdf))  
  df <- rbind(df, subdf)
  row_id_counter <- max(df$row_id, na.rm = TRUE) + 1 
}

# Compute geometric mean 
#geom_mean <- data.frame(lambda = as.numeric(colnames(df)[1:(ncol(df) - 2)]),
#                        value = colMeans(-log10(df[, 1:(ncol(df) - 2)]), na.rm = TRUE))
#geom_mean$value[1] <- NA
#geom_mean$value[nrow(geom_mean)] <- NA
#geom_mean$context <- "geometric mean"


df$PRS <- rep(dimnames(pval_array)[[1]], length(envs))

df_long <- melt(df, 
                id.vars = c("context", "row_id", "PRS"), 
                variable.name = "lambda", 
                value.name = "value")

df_long$lambda <- as.numeric(as.character(df_long$lambda))
df_long$context <- as.factor(df_long$context)
levels(df_long$context)[levels(df_long$context) == "sex"] <- "Sex"
levels(df_long$context)[levels(df_long$context) == "age"] <- "Age"
levels(df_long$context)[levels(df_long$context) == "Alcohol_intake_frequency"] <- "Alcohol intake frequency"
levels(df_long$context)[levels(df_long$context) == "Smoking_status"] <- "Smoking status"

top_prs <- df_long %>%
  group_by(context, PRS) %>%
  summarise(max_val = max(-log10(value), na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(order_by = max_val, n = 1, with_ties = FALSE)

# Keep up to 3 PGS, only if they are significant
labels <- df_long %>%
  group_by(context, PRS) %>%
  summarise(max_val = max(-log10(value), na.rm = TRUE), 
            lambda_at_max = lambda[which.max(-log10(value))],
            .groups = "drop") %>%
  filter(max_val > -log10(threshold)) %>%                     
  group_by(context) %>%
  slice_max(order_by = max_val, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  select(context, PRS, lambda = lambda_at_max, value = max_val)

format_pgs_label <- function(s) {
  s <- gsub("_", " ", s)
  s <- capitalize_first(s)
  paste0(s, " PGS")
}

labels$PRS <- sapply(labels$PRS, format_pgs_label)
df_long[df_long$PRS=="height",]
if (length(envs) > 1) {
png(paste0("/Users/manuelacostantino/Desktop/blue_lines",pheno,".png"), width = 8, height = 5, res=200, units = "in")
ggplot(df_long, aes(x = lambda, y = -log10(value), group = row_id)) +
  geom_line(color = "#00CCFF") +
  geom_hline(yintercept = -log10(threshold), 
             linetype = "dashed", color = "black", size = 0.6) +
  geom_text_repel(
    data = labels,
    aes(x = lambda, y = value, label = PRS),
    size = 2, color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    #title = paste("PRSxC effects for", gsub("_", " ", pheno)),
    x = "Lambda",
    y = expression(-log[10](p))
  ) +
  facet_wrap(~ context, scales = "free_y", ncol = 3) +
  #coord_cartesian(ylim = c(0, max(-log10(df_long$value))+2)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "white")
  )
dev.off()
} else {
  png(paste0("/Users/manuelacostantino/Desktop/blue_lines_one",pheno,".png"), width = 3, height = 4, res=200, units = "in")
  ggplot(df_long, aes(x = lambda, y = -log10(value), group = row_id)) +
    geom_line(color = "#00CCFF") +
    geom_hline(yintercept = -log10(threshold), 
               linetype = "dashed", color = "black", size = 0.6) +
    geom_text_repel(
      data = labels,
      aes(x = lambda, y = value, label = PRS),
      size = 3, color = "black",
      inherit.aes = FALSE
    ) +
    labs(
      #title = paste("PRSxC effects for", gsub("_", " ", pheno)),
      x = "Lambda",
      y = expression(-log[10](p))
    ) +
    #coord_cartesian(ylim = c(0, max(-log10(df_long$value))+2)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.2),
      panel.background = element_rect(fill = "white")
    )
  dev.off()
}
