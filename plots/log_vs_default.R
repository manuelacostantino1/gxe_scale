library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)

# Load in PRS for testosterone and height, the two phenos and their log versions, sex and age
height <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Height/Height674178.pheno", header=TRUE)
testosterone <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Testosterone/Testosterone674178.pheno", header=TRUE)
covars <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE)
#prs_height <- read.table("/gpfs/data/ukb-share/dahl/manuela/ukb_prs_downloading/polygenic_scores/sPRS_height/sPRS_height678224.pheno", header=TRUE)
#prs_bmi <- read.table("/gpfs/data/ukb-share/dahl/manuela/ukb_prs_downloading/polygenic_scores/sPRS_BMI/sPRS_BMI678224.pheno", header=TRUE)

# Get best thresholds for both height and bmi
tes_tfile <- read.table("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out/processed_testosterone_whitebrit_PGS.txt")
height_tfile <- read.table("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out/processed_height_whitebrit_PGS.txt")
tes_thresh <- format(tes_tfile[1], scientific=FALSE)
height_thresh <- format(height_tfile[1], scientific=FALSE)

# Open prs files
prs_tes <- read.table(paste("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/testosterone/testosterone_pgs.",tes_thresh,".sscore", sep=""), header=TRUE, strip.white = TRUE, comment.char = "")
colnames(prs_tes) <- c("FID", "IID", "allele_tes", "dosage_tes", "avg_score_tes", "sum_score_tes")
prs_height <- read.table(paste("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/height/height_pgs.",height_thresh,".sscore", sep=""), header=TRUE, strip.white = TRUE, comment.char = "")
colnames(prs_height) <- c("FID", "IID", "allele_height", "dosage_height", "avg_score_height", "sum_score_height")

merged_df <- height %>%
  merge(testosterone, by = "IID") %>%
  merge(covars, by = "IID") %>%
  merge(prs_height, by = "IID") %>%
  merge(prs_tes, by = "IID")

#head(merged_df)

df <- data.frame(height=merged_df$X50, testosterone=merged_df$X30850, log_height=log(merged_df$X50), log_testosterone=log(merged_df$X30850), sex=as.factor(merged_df$X31.0.0), age=merged_df$X21003.0.0, prs_height=merged_df$sum_score_height, prs_tes=merged_df$sum_score_tes)
df <- na.omit(df)
levels(df$sex) <- c("Female", "Male")
head(df)

df_long <- df %>%
  pivot_longer(
    cols = c(height, log_height, testosterone, log_testosterone),
    names_to = "variable",
    values_to = "value"
  )

# Remove outliers per variable (1st and 99th percentile)
df_long_filtered <- df_long %>%
  group_by(variable) %>%
  mutate(
    lower = quantile(value, 0.001, na.rm = TRUE),
    upper = quantile(value, 0.999, na.rm = TRUE),
    value = ifelse(value < lower | value > upper, NA, value)
  ) %>%
  ungroup() %>%
  filter(!is.na(value))

# Rename variables for nicer facet labels
df_long_filtered <- df_long_filtered %>%
  mutate(variable = recode(variable,
                           height = "Height",
                           log_height = "log(Height)",
                           testosterone = "Testosterone",
                           log_testosterone = "log(Testosterone)"))

binwidth_value <- 1

png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/test_vs_log_test_hist_4panels.png", width = 1300, height = 1200)

ggplot(df_long_filtered, aes(x = value, fill = sex)) +
  geom_density(alpha = 0.8) +
  facet_wrap(~ variable, scales = "free", ncol = 2) +  # 1 col for vertical panels
  labs(x = "Value", y = "Density", fill = "Sex") +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 33),
    axis.text = element_text(size = 26),
    strip.text = element_text(size = 30),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 30),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = c("Male" = "#66CCFF", "Female" = "#660099"))

dev.off()

# Filter for testosterone-related variables
df_testosterone <- df_long_filtered %>%
  filter(variable %in% c("Testosterone", "log(Testosterone)"))

# Filter for height-related variables
df_height <- df_long_filtered %>%
  filter(variable %in% c("Height", "log(Height)"))

# Plot for testosterone variables
#png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/testosterone_hist_horizontal.png", width = 1300, height = 600)
png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/testosterone_hist_vertical.png",width = 3, height = 4, units = "in", res= 300)
ggplot(df_testosterone, aes(x = value, fill = sex)) +
  geom_density(alpha = 0.8, size = 0.3) +  # thinner lines for smaller figure
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  labs(x = "Value", y = "Density", fill = "Sex") +
  theme_minimal(base_size = 8) +  # smaller base font
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Male" = "#66CCFF", "Female" = "#660099"))
dev.off()

# Plot for height variables
#png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/height_hist_horizontal.png", width = 1300, height = 600)
png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/height_hist_vertical.png",width = 3, height = 4, units = "in", res= 300)
ggplot(df_height, aes(x = value, fill = sex)) +
  geom_density(alpha = 0.8, size = 0.3) +  # thinner lines for smaller figure
  facet_wrap(~ variable, scales = "free", nrow = 2) +
  labs(x = "Value", y = "Density", fill = "Sex") +
  theme_minimal(base_size = 8) +  # smaller base font
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Male" = "#66CCFF", "Female" = "#660099"))
dev.off()


#####################

# this part makes log vs default plots for height and testosterone

# Make sure no log problems
df <- df %>%
  mutate(
    log_height = log(height),
    log_testosterone = log(testosterone)
  )

# Sample 50,000 points max to avoid crashing
set.seed(1)
df_plot <- df %>%
  filter(
    !is.na(height), !is.na(log_height),
    !is.na(testosterone), !is.na(log_testosterone),
    is.finite(log_height), is.finite(log_testosterone)
  ) %>%
  sample_n(min(50000, n()))

p_height <- ggplot(df_plot, aes(x = height, y = log_height, color = sex)) +
  geom_point(alpha = 0.6, size = 3) +  # shape 16 is default solid circle
  scale_color_manual(values = c("Male" = "#66CCFF", "Female" = "#660099")) +
  labs(x = "Height", y = "log(Height)", color = "Sex") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

p_tes <- ggplot(df_plot, aes(x = testosterone, y = log_testosterone, color = sex)) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_manual(values = c("Male" = "#66CCFF", "Female" = "#660099")) +
  labs(x = "Testosterone", y = "log(Testosterone)", color = "Sex") +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )


# Arrange side-by-side
png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/log_vs_def.png", width = 13, height = 10, units = "in", res= 300)
grid.arrange(p_tes, p_height, ncol = 2)
dev.off()
png(
  "/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/height_vs_log_height.png",
  width = 3, height = 3, units = "in", res = 300
)
ggplot(df_plot, aes(x = height, y = log_height, color = sex)) +
  geom_point(alpha = 0.6, size = 1.5) +  # smaller points for smaller canvas
  scale_color_manual(values = c("Male" = "#66CCFF", "Female" = "#660099")) +
  labs(x = "Height", y = "log(Height)", color = "Sex") +
  theme_minimal(base_size = 8) +  # smaller base font
  theme(
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3)
  )
dev.off()

# Testosterone plot
png(
  "/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/testosterone_vs_log_testosterone.png",
  width = 3, height = 3, units = "in", res = 300
)
ggplot(df_plot, aes(x = testosterone, y = log_testosterone, color = sex)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Male" = "#66CCFF", "Female" = "#660099")) +
  labs(x = "Testosterone", y = "log(Testosterone)", color = "Sex") +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.3)
  )
dev.off()


