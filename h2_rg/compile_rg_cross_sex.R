all_phenos <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Arm_fat-free_mass_left", "Creatinine", "Triglycerides", 
    "Arm_fat-free_mass_right", "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "Basophill_count", "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol", "Eosinophill_count", "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin")

parse_rg <- function(filepath){
    # read the file
    lines <- readLines(filepath)

    # get the rg and standard error
    rg_line <- lines[length(lines) - 3]
    df <- read.table(text = rg_line, header = FALSE, stringsAsFactors = FALSE)
    rg <- df$V3
    rg_se <- df$V4

    # get the two heritabilities
    h2_female_line <- lines[length(lines) - 24]
    nums <- regmatches(h2_female_line, gregexpr("[0-9\\.]+", h2_female_line))[[1]]
    h2_female <- as.numeric(nums[2])
    h2_female_se  <- as.numeric(nums[3])

    h2_male_line <- lines[length(lines) - 32]
    nums <- regmatches(h2_male_line, gregexpr("[0-9\\.]+", h2_male_line))[[1]]
    h2_male <- as.numeric(nums[2])
    h2_male_se  <- as.numeric(nums[3])

    # return the values
    return(list(rg = rg, rg_se = rg_se,
                h2_female = h2_female, h2_female_se = h2_female_se,
                h2_male = h2_male, h2_male_se = h2_male_se))
}


results <- data.frame()

for (pheno in all_phenos) {
    log_path <- paste0("/gpfs/data/ukb-share/dahl/manuela/ldsc_results/cross_sex/", tolower(pheno), ".log")
    log_log_path <- paste0("/gpfs/data/ukb-share/dahl/manuela/ldsc_results/cross_sex_log/", tolower(pheno), ".log")

    if (file.exists(log_path) & file.exists(log_log_path)) {
        parsed_default <- parse_rg(log_path)
        parsed_log <- parse_rg(log_log_path)

        results <- rbind(
            results,
            data.frame(
                phenotype = pheno,

                # default
                rg_default = parsed_default$rg,
                rg_se_default = parsed_default$rg_se,
                h2_female_default = parsed_default$h2_female,
                h2_female_se_default = parsed_default$h2_female_se,
                h2_male_default = parsed_default$h2_male,
                h2_male_se_default = parsed_default$h2_male_se,

                # log
                rg_log = parsed_log$rg,
                rg_se_log = parsed_log$rg_se,
                h2_female_log = parsed_log$h2_female,
                h2_female_se_log = parsed_log$h2_female_se,
                h2_male_log = parsed_log$h2_male,
                h2_male_se_log = parsed_log$h2_male_se,

                stringsAsFactors = FALSE
            )
        )
    } else {
        message("Missing file(s) for ", pheno)
    }
}

write.csv(results, file = "/gpfs/data/ukb-share/dahl/manuela/scale_project/rgs_cross_sex.csv", quote = FALSE, row.names = FALSE)