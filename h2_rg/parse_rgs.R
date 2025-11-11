
all_phenos <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Arm_fat-free_mass_left", "Creatinine", "Triglycerides", 
    "Arm_fat-free_mass_right", "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "Basophill_count", "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol", "Eosinophill_count", "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin")

# Function to parse rg output
parse_rg <- function(filepath){

    # read the file
    lines <- readLines(filepath)

    # get the rg and standard error
    rg_line <- lines[length(lines) - 3]
    df <- read.table(text = rg_line, header = FALSE, stringsAsFactors = FALSE)
    rg <- df$V3
    rg_se <- df$V4

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
    return(list(rg = rg, rg_se = rg_se,
                      h2_default = h2_default, h2_default_se = h2_default_se,
                      h2_log = h2_log, h2_log_se = h2_log_se))
}


results <- data.frame(
  pheno = character(),
  rg = numeric(),
  rg_se = numeric(),
  h2_default = numeric(),
  h2_default_se = numeric(),
  h2_log = numeric(),
  h2_log_se = numeric(),
  stringsAsFactors = FALSE
)

for (pheno in all_phenos) {
    filepath <- paste0("/gpfs/data/ukb-share/dahl/manuela/ldsc_results/cross_scale/", tolower(pheno), ".log")
    
    if (file.exists(filepath)) {
        vals <- parse_rg(filepath)
        
        results <- rbind(
          results,
          data.frame(
            pheno = pheno,
            rg = vals$rg,
            rg_se = vals$rg_se,
            h2_default = vals$h2_default,
            h2_default_se = vals$h2_default_se,
            h2_log = vals$h2_log,
            h2_log_se = vals$h2_log_se,
            stringsAsFactors = FALSE
          )
        )
    } else {
        warning(paste("File not found:", filepath))
    }
}

write.table(results, file = "/gpfs/data/ukb-share/dahl/manuela/ldsc_results/rg_results_cross_scale.txt", sep = "\t", row.names = FALSE, quote = FALSE)