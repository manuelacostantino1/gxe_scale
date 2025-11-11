library(VennDiagram)
library(ggplot2)
library(tidyr)
library(dplyr)

df <- read.table("/gpfs/data/ukb-share/dahl/manuela/scale_project/gwas_sharing.txt", stringsAsFactors=FALSE) 
colnames(df) <- c("pheno", "default_union", "log_union", "all_three")
head(df)

# Venn diagram for testosterone and height
for (pheno in c("testosterone", "height")){
png(paste0("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/even_venn/descaling_gwas_clumped_venn_", pheno, ".png"),
    width = 800, height = 800)
    inverted <- (as.numeric(df[df$pheno == pheno, "default_union"])  < as.numeric(df[df$pheno == pheno, "log_union"]))

if (!inverted) {
    area1_use <- df[df$pheno == pheno, "default_union"] 
    area2_use <- df[df$pheno == pheno, "log_union"]
    cross_use <- df[df$pheno == pheno, "all_three"]
    categories <- c("Default scale", "Log scale")
} else {
    # Swap areas
    area1_use <- df[df$pheno == pheno, "log_union"]
    area2_use <- df[df$pheno == pheno, "default_union"] 
    cross_use <- df[df$pheno == pheno, "all_three"]
    categories = c("Log scale", "Default scale")
}
grid.newpage()
venn.plot <- draw.pairwise.venn(
    area1 = area1_use,
    area2 = area2_use,
    cross.area = cross_use,
    category = categories,
    fill = c("lightsteelblue1", "navyblue"),  
    col  = c("steelblue3", "midnightblue"),
    alpha = 0.5,
    lty = 1,
    lwd = 2,
    cex = 4,
    fontfamily = "sans",
    cat.cex = 3,
    cat.fontfamily = "sans",
    cat.pos = c(0, 0),
    cat.dist = c(0.07, 0.07),
    ind = FALSE,
    scaled = FALSE
)
pushViewport(viewport(y = 0.5, height = 1, width = 1))  
grid.draw(gTree(children = venn.plot))
popViewport()
dev.off()
}

# Calculate shared vs specific
df$specific_default <- as.numeric(df$default_union) - as.numeric(df$all_three)
df$specific_log <- as.numeric(df$log_union) - as.numeric(df$all_three)
df$pheno <- paste0(toupper(substring(df$pheno, 1, 1)), substring(df$pheno, 2)) # Capitalize first letter
df$pheno <- gsub("_", " ", df$pheno) # Replace underscores with spaces
df <- df[!(df$pheno == "Whradjbmi emdin"), ]
df <- df[!(df$pheno == "Arm fat-free mass left"), ]
df <- df[!(df$pheno == "Arm fat-free mass right"), ]
df <- df[!(df$pheno == "Ea4"), ]
df$pheno[df$pheno == "Arm fat-free mass avg"] <- "Arm fat-free mass"
df$pheno[df$pheno == "Whradjbmi zhu"] <- "WHRadjBMI"
df$pheno[df$pheno == "Systolicbp auto"] <- "Systolic BP"
df$pheno[df$pheno == "Diastolicbp auto"] <- "Diastolic BP"

# Reformat for ggplot
df_long <- df %>%
  select(pheno, specific_default, specific_log, all_three) %>%
  pivot_longer(cols = c(specific_default, specific_log, all_three),
               names_to = "category", values_to = "count")

df_long$category <- factor(df_long$category,
                           levels = c("specific_default", "all_three", "specific_log"),
                           labels = c("Specific to default", "Shared", "Specific to log"))

df$total_hits <- df$specific_default + df$specific_log + df$all_three
df_long$pheno <- factor(df_long$pheno,
                        levels = df$pheno[rev(order(df$total_hits, decreasing = TRUE))])

png(paste0("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/even_venn/all_phenos_relative_gwas.png"),
    width = 8, height = 5, unit="in" ,  res = 600)
ggplot(df_long, aes(x = pheno, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Phenotype",
    y = "Number of significant loci",
    fill = "Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Specific to default" = "darkolivegreen4",
    "Shared" = "darkseagreen3",
    "Specific to log" = "darkslategray4"
  ))
dev.off()