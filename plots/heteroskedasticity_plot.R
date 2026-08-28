library(ggplot2)
library(ggrepel)
library(tidyr)
library(reshape)
library(dplyr)

# List of phenos
all_phenos <- c(
  "Birth_weight", "HDL", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
  "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Creatinine", "Triglycerides",
  "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
  "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "LDL",  "SHBG", "BMI", "Cholesterol", "Hip_circumference",
  "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu"
)

# "Leukocyte_count" is missing

quantnorm <- function(Y) {
  obs <- which(!is.na(Y))
  ranks <- rank(Y[obs])
  quantiles <- qnorm((ranks + 1) / (length(ranks) + 2))
  if (!all.equal(ranks, rank(quantiles))) stop("Bad rank problem")
  Y[obs] <- quantiles
  scale(Y)
}


# Open covars and get ID and sex
system("dx download /phenotypes/extracted_phenotypes/covariates/covars_age2.pheno")
covars <- read.table("covars_age2.pheno", header=TRUE)
pheno_df <- covars[,1:3]

# Open all phenos and make a df
for (pheno in all_phenos){
  local_file <- paste0("phenos/", pheno, ".pheno")
  if (!file.exists(local_file)) {
    pheno_path <- paste0("/phenotypes/extracted_phenotypes/", pheno, "/", pheno, ".pheno")
    system(paste0("dx download ", pheno_path, " -o ", local_file))
  }
  pheno_dat <- read.table(local_file, header=TRUE)
  colnames(pheno_dat)[3] <- pheno
  pheno_df <- merge(pheno_df, pheno_dat, by=c("FID","IID"))
}

head(pheno_df)

# Calculate heteroskedasticity for the default scale, log scale and RINT scale
results <- data.frame()

for (pheno in all_phenos) {
  
  y <- pheno_df[[pheno]]
  sex <- pheno_df$X31.0.0.x   # 1 = male, 0 = female
  
  # Default
  r_default <- var(y[sex == 1], na.rm = TRUE) / var(y[sex == 0], na.rm = TRUE)
  
  # Log (only if valid)
  r_log <- if (all(y > 0, na.rm = TRUE)) {
    var(log(y)[sex == 1], na.rm = TRUE) / var(log(y)[sex == 0], na.rm = TRUE)
  } else {
    NA
  }
  
  # Quantile normalized
  y_qn <- as.numeric(quantnorm(y))
  r_qn <- var(y_qn[sex == 1], na.rm = TRUE) / var(y_qn[sex == 0], na.rm = TRUE)
  
  results <- rbind(results, data.frame(
    phenotype = pheno,
    default_ratio = r_default,
    log_ratio = r_log,
    quantnorm_ratio = r_qn
  ))
}

head(results)

# Make plot and save
results <- results[results$phenotype != "Testosterone", ]
plot_df <- results |>
  pivot_longer(
    cols = c(log_ratio, quantnorm_ratio),
    names_to = "scale",
    values_to = "ratio"
  )
plot_df$scale <- recode(plot_df$scale,
                        log_ratio = "Log scale",
                        quantnorm_ratio = "Quantile normalized"
)



p <- ggplot(plot_df, aes(x = default_ratio, y = ratio, label = phenotype)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_text_repel(
    size = 3,
    max.overlaps = 5
  ) +
  facet_wrap(~ scale, nrow = 1) +
  theme_classic() +
  labs(
    x = "Male / Female variance ratio (default scale)",
    y = "Male / Female variance ratio (transformed scale)"
  )
ggsave("plot.png", p, width = 10, height = 4, dpi = 600)
