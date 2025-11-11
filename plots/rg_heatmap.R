library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggnewscale)

LDSC_RESULTS <- "/gpfs/data/ukb-share/dahl/manuela/ldsc_results/"
pheno <- "height"

# Write parsing function
parse_rg <- function(pheno, method1, method2){

    # read the file
    filepath <- paste0(LDSC_RESULTS, pheno, "/", method1, "_", method2, "_", pheno, ".log")
    lines <- readLines(filepath)

    # get the rg and standard error
    # LDSC outputs genetic correlation in a table near the end, often 4th line from last
    rg_line <- lines[length(lines) - 3]
    df <- read.table(text = rg_line, header = FALSE, stringsAsFactors = FALSE)
    rg <- as.numeric(df$V3)
    rg_se <- as.numeric(df$V4)

    # get the two heritabilities
    h2_default_line <- lines[length(lines) - 24]
    nums <- regmatches(h2_default_line, gregexpr("[0-9\\.]+", h2_default_line))[[1]]
    h2_default <- as.numeric(nums[2])
    h2_default_se  <- as.numeric(nums[3])

    h2_log_line <- lines[length(lines) - 32]
    nums <- regmatches(h2_log_line, gregexpr("[0-9\\.]+", h2_log_line))[[1]]
    h2_log <- as.numeric(nums[2])
    h2_log_se  <- as.numeric(nums[3])

    # return the values
    return(list(method1 = as.factor(method1), method2 = as.factor(method2),
                rg = rg, rg_se = rg_se,
                h2_default = h2_default, h2_default_se = h2_default_se,
                h2_log = h2_log, h2_log_se = h2_log_se))
}

# Read in results
methods <- c("default", "log", "logmales", "males",  "females", "logfemales")
results <- list()  # use a list first for efficiency

for (method1 in methods){
    for (method2 in methods){
        if (method1 == method2) next
        print(paste0("Parsing ", method1, " and ", method2))
        res <- parse_rg(pheno, method1, method2)
        results[[paste(method1, method2, sep = "_")]] <- res
    }
}

# Convert list to data frame
results_df <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))

results_df
results_df$method1 <- as.character(results_df$method1)
results_df$method2 <- as.character(results_df$method2)
n <- length(methods)

# Initialize matrices
h2_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(methods, methods))
rg_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(methods, methods))

# Get the rgs
for(i in 1:nrow(results_df)){
  m1 <- results_df$method1[i]
  m2 <- results_df$method2[i]
  
  # Fill lower triangle only
  rg_matrix[m2, m1] <- results_df$rg[i]
}

# Get heritabilities
for(m in methods){
  idx <- which(results_df$method1 == m | results_df$method2 == m)[1]
  h2_matrix[m, m] <- results_df$h2_default[idx]
}

h2_long <- melt(h2_matrix, varnames = c("Trait1","Trait2"), value.name = "h2") %>%
  filter(!is.na(h2))
rg_long <- melt(rg_matrix, varnames = c("Trait1","Trait2"), value.name = "rg") %>%
  filter(!is.na(rg)) %>%
  filter(match(Trait1, methods) > match(Trait2, methods))


h2_long$Trait1 <- factor(h2_long$Trait1, levels = methods)
h2_long$Trait2 <- factor(h2_long$Trait2, levels = methods)
rg_long$Trait1 <- factor(rg_long$Trait1, levels = methods)
rg_long$Trait2 <- factor(rg_long$Trait2, levels = methods)

method_labels <- c(
  "default"    = "All - Default",  
  "log"        = "All - Log",
  "males"      = "Males - Default",
  "logmales"   = "Males - Log",
  "females"    = "Females - Default",
  "logfemales" = "Females - Log"
)


rg_long$Trait1 <- factor(rg_long$Trait1, levels = methods, labels = method_labels)
rg_long$Trait2 <- factor(rg_long$Trait2, levels = methods, labels = method_labels)
h2_long$Trait1 <- factor(h2_long$Trait1, levels = methods, labels = method_labels)
h2_long$Trait2 <- factor(h2_long$Trait2, levels = methods, labels = method_labels)
rg_long <- rg_long %>% filter(as.numeric(Trait1) > as.numeric(Trait2))
h2_long <- h2_long %>% filter(Trait1 == Trait2)


png(paste0("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/rg_heatmap_", pheno, ".png"),
    width = 3, height = 4, units = "in", res = 600)

ggplot() +
  # Genetic correlations
  geom_tile(data = rg_long, aes(x = Trait2, y = Trait1, fill = rg)) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "darkslateblue",
    midpoint = 0, limits = c(-0.05, 1.02),
    name = "Genetic correlation",
    guide = guide_colorbar(
      label = TRUE, barwidth = 5.2, barheight = 1,
      ticks = TRUE, ticks.linewidth = 0.5,
      nbin = 100, title.position = "top", title.hjust = 0.5
    )
  ) +
  # New fill for heritabilities
  ggnewscale::new_scale_fill() +
  geom_tile(data = h2_long, aes(x = Trait2, y = Trait1, fill = h2)) +
  scale_fill_gradient(
    low = "gray79", high = "gray23",
    limits = c(0, 0.5), name = "Heritability",
    guide = guide_colorbar(
      label = TRUE, barwidth = 5.2, barheight = 1,
      ticks = TRUE, ticks.linewidth = 0.5,
      title.position = "top", title.hjust = 0.5
    )
  ) +
  # Add text labels
  geom_text(data = rg_long, aes(x = Trait2, y = Trait1, label = round(rg, 3)),
            color = "white", size = 2.5, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = h2_long, aes(x = Trait2, y = Trait1, label = round(h2, 3)),
            color = "white", size = 2.5, fontface = "bold", inherit.aes = FALSE) +
  # Fix y-axis order
  scale_y_discrete(limits = rev(levels(rg_long$Trait1))) +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold", color = "black"),
    axis.text.y = element_text(size = 7, face = "bold", color = "black"),
    axis.title = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.box.just = "left",
    legend.title.position = "top"
  )

dev.off()