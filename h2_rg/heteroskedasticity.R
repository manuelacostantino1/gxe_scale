all_phenos <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Arm_fat-free_mass_left", "Creatinine", "Triglycerides", 
    "Arm_fat-free_mass_right", "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "Basophill_count", "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol", "Eosinophill_count", "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin")

# Open covar file
covars <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE)

# Start output df
df <- data.frame(pheno = character(), var_female = numeric(), var_male = numeric(), stringsAsFactors = FALSE)

for (pheno in all_phenos) {
    
    # Open  pheno file
    pheno_path <- paste("/gpfs/data/ukb-share/extracted_phenotypes/", pheno, sep="")
    pheno_file <- list.files(pheno_path, full.names = TRUE)[1]
    pheno_table <- read.table(pheno_file, header = TRUE)
    colnames(pheno_table) <- c("FID", "IID", "pheno")    

    # Merge with covar file
    merged_df <- merge(pheno_table, covars, by = c("IID", "FID"))
    merged_df$log_pheno <- log(merged_df$pheno)

    # Append sex-specific variances to output df
    males <- merged_df[merged_df$X31.0.0==1,]
    females <- merged_df[merged_df$X31.0.0==0,]
    
    var_fem_def <- var(females$pheno, na.rm=TRUE)
    var_mal_def <- var(males$pheno, na.rm=TRUE)
    var_fem_log <- var(females$log_pheno, na.rm=TRUE)
    var_mal_log <- var(males$log_pheno, na.rm=TRUE)

    new_row <- data.frame(pheno = pheno, var_female = var_fem_def, var_male   = var_mal_def, log_var_female = var_fem_log,log_var_male   = var_mal_log, stringsAsFactors = FALSE)

    df <- rbind(df, new_row)
}

# Write output df to file
head(df)
write.table(df, file="/gpfs/data/ukb-share/dahl/manuela/scale_project/heteroskedasticity_all_phenos.txt", sep="\t", quote=FALSE, row.names=FALSE)