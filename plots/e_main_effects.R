library(dplyr)
library(tidyr)
library(ggplot2)

# Make a list of the phenos and the envs and matrix to save the effect sizes
all_phenos <- c("Birth_weight", "HDL", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Creatinine", "Triglycerides", 
     "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D",  "LDL", "Leukocyte_count", "SHBG", "BMI", "Cholesterol",  "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu" )
envs <- c("sex", "age", "Smoking_status", "Alcohol_intake_frequency", "Statins")

# Data frame to save effect sizes
results <- data.frame(pheno=character(), beta_sex=numeric(), beta_age=numeric(), beta_smoking1=numeric(),beta_smoking1=numeric(), beta_alcohol=numeric(), beta_statins=numeric(), stringsAsFactors=FALSE)

# Open all the envs and make a df with all of them
cov_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE)
df <- data.frame(FID = cov_table$FID, IID = cov_table$IID, sex = as.factor(cov_table$X31.0.0), age = scale(as.numeric(cov_table$X21003.0.0)), center = as.factor(cov_table$X54.0.0))
statins_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Statins/Statins674178.pheno", header=TRUE)
colnames(statins_table)[3] <- "statins"
smoking_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Smoking_status/Smoking_status674178.pheno", header=TRUE)
colnames(smoking_table)[3] <- "smoking"
alcohol_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/Alcohol_intake_frequency/Alcohol_intake_frequency674178.pheno", header=TRUE)
colnames(alcohol_table)[3] <- "alcohol"
df <- merge(df, statins_table[,c("FID", "IID", "statins")], by=c("FID", "IID"), all.x=TRUE)
df <- merge(df, smoking_table[,c("FID", "IID", "smoking")], by=c("FID", "IID"), all.x=TRUE)
df <- merge(df, alcohol_table[,c("FID", "IID", "alcohol")], by=c("FID", "IID"), all.x=TRUE)
df$statins <- as.factor(df$statins)
df$smoking <- as.factor(df$smoking)
df$alcohol <- scale(as.numeric(df$alcohol))

# Iterate through phenos and for each make a df
for (my_pheno in all_phenos) {

    # Iterate through envs and for each get the effect size and save it to the matrix
    pheno_path <- paste("/gpfs/data/ukb-share/extracted_phenotypes/", my_pheno, sep="")
    pheno_file <- list.files(pheno_path, full.names = TRUE)[1]
    pheno_table <- read.table(pheno_file, header = TRUE)
    df_temp <- merge(df, pheno_table, by=c("FID", "IID"), all.x=TRUE)
    colnames(df_temp)[which(colnames(df_temp) == colnames(pheno_table)[3])] <- "pheno"
    print(my_pheno)

    mod <- lm(scale(pheno) ~ sex + age + smoking + alcohol + statins + center, data=df_temp)
    beta_sex <- summary(mod)$coefficients["sex1", "Estimate"]
    beta_age <- summary(mod)$coefficients["age", "Estimate"]
    beta_smoking1 <- summary(mod)$coefficients["smoking1", "Estimate"]
    beta_smoking2 <- summary(mod)$coefficients["smoking2", "Estimate"]
    beta_alcohol <- summary(mod)$coefficients["alcohol", "Estimate"]
    beta_statins <- summary(mod)$coefficients["statins1", "Estimate"]
    if (my_pheno =="SystolicBP_auto") my_pheno <- "Systolic BP"
    if (my_pheno =="DiastolicBP_auto") my_pheno <- "Diastolic BP"
    if (my_pheno =="Arm_fat-free_mass_avg") my_pheno <- "Arm_fat-free_mass"
    if (my_pheno =="WHRadjBMI_Zhu") my_pheno <- "WHRadjBMI"
    my_pheno <- gsub("_", " ", my_pheno)

    results <- rbind(results, data.frame(pheno = my_pheno, beta_sex = beta_sex, beta_age = beta_age, beta_smoking1 = beta_smoking1,beta_smoking2 = beta_smoking2, beta_alcohol = beta_alcohol, beta_statins = beta_statins))
}
print(results)

# Make heatmap plot
results_long <- results %>%
  pivot_longer(
    cols = starts_with("beta_"),
    names_to = "environment",
    values_to = "effect"
  ) %>%
  mutate(environment = gsub("beta_", "", environment)) %>%
  mutate(environment = dplyr::recode(environment,
      sex = "Sex",
      age = "Age",
      smoking1 = "Smoking: Former",
      smoking2 = "Smoking: Current",
      alcohol = "Alcohol Intake Frequency",
      statins = "Statin Use"
  ))

results_long$pheno <- factor(results_long$pheno,levels = rev(unique(results_long$pheno)))

png("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/e_main_effects_heatmap.png", width = 3000, height = 2400, res=600)
ggplot(results_long, aes(x = environment, y = pheno, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 8) +   
  theme(
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.title = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Environment",
    y = "Phenotype",
    title = "Effect of Environments on Phenotypes"
  )
dev.off()