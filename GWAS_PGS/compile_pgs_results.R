rm(list=ls())

home_dir <- "/Users/manuelacostantino/Documents/descaler_inputs/"
pheno_list <- c("Alcohol_intake_frequency674178", "Birth_weight674178", "BMI674178", "Calcium674178", "Cholesterol674178", "Creatinine674178",  "Glucose674178", "HDL674178", "HbA1c674178", "Height674178", "Hip_circumference674178", "IGF-1674178", "LDL674178", "Leukocyte_count674178", "Mean_corpuscular_volume674178", "Platelet_count674178", "Pulse_rate674178", "RBC674178", "SHBG674178", "SystolicBP_auto674178", "Testosterone674178", "Triglycerides674178", "Urate674178", "Urea674178", "Vitamin_D674178", "Waist_circumference674178", "Whole_body_fat_mass674178", "WHRadjBMI_Zhu",  "EA4", "C-reactive_protein674178", "Hip_to_waist674178", "Arm_fat-free_mass_avg674178")
#"DiastolicBP_auto674178", "FEV1674178",
valid_pops <- c("white_euro")
# Diastolic bp , ea4 and fev1 crsaahing
tau <- c(0.0000000001, 0.00000001, 0.000001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5)

st <- format(Sys.time(), "%Y-%m-%d")
base_cols <- c(  '#6699CC', '#EE99AA','#994455','#004488','#DDDDDD' )
trans_cols <- c(paste0(base_cols[1],"99"),paste0(base_cols[2],"99"),paste0(base_cols[3],"99"),
                paste0(base_cols[4],"99"),paste0(base_cols[5],"99"))

output_all <- array(NA, dim=c(length(pheno_list),  7, length(valid_pops)),dimnames=list(pheno_list, c("thresh0", "thresh1", "r20", "r21", "ci_95_l","ci_95_u", "pval"), valid_pops))
output_fem <- array(NA, dim=c(length(pheno_list),  7, length(valid_pops)),dimnames=list(pheno_list, c("thresh0", "thresh1", "r20", "r21", "ci_95_l", "ci_95_u", "pval"), valid_pops))
output_male <- array(NA, dim=c(length(pheno_list), 7, length(valid_pops)), dimnames=list(pheno_list,  c("thresh0", "thresh1", "r20", "r21", "ci_95_l", "ci_95_u", "pval"), valid_pops))


for(pheno in pheno_list){
  if (grepl("FEV1674178", pheno)) {
    phenoNoDigits <- "FEV1"
  } else if (grepl("IGF-1674178", pheno)) {
    phenoNoDigits <- "IGF-1"
  } else if (grepl("EA4", pheno)) {
    phenoNoDigits <- "EA4"
  } else {
    phenoNoDigits <- gsub("[0-9]+$", "", pheno)
  }
  phenoLower <- tolower(phenoNoDigits)
  
  pgs_white_euro_all <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_white_euro_all.txt"), header=TRUE)
  #pgs_afr_all <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_afr_all.txt"), header=TRUE)
  #pgs_asn_all <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_asn_all.txt"), header=TRUE)
  pgs_white_euro_fem <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_white_euro_female.txt"), header=TRUE)
  #pgs_afr_fem <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_afr_female.txt"), header=TRUE)
  #pgs_asn_fem <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_asn_female.txt"), header=TRUE)
  pgs_white_euro_male <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_white_euro_male.txt"), header=TRUE)
  #pgs_afr_male <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_afr_male.txt"), header=TRUE)
  #pgs_asn_male <- read.table(paste0(home_dir,"log_vs_lin/log_vs_lin_",phenoLower,"_asn_male.txt"), header=TRUE)
  
  
  add_line <- function(pheno, pop ,output, table){
    output[pheno,"thresh0" ,pop] <- table[1,"thresh0"]
    output[pheno,"thresh1" ,pop] <- table[1,"thresh1"]
    output[pheno,"r20" ,pop] <- table[1,"r20"]
    output[pheno,"r21" ,pop] <- table[1,"r21"]
    output[pheno,"ci_95_l" ,pop] <- table[1,"ci_95"]
    output[pheno,"ci_95_u" ,pop] <- table[2,"ci_95"]
    output[pheno,"pval" ,pop] <- table[1,"pval"]
    return(output)
  }
  
  output_all <-add_line(pheno, "white_euro",  output_all, pgs_white_euro_all)
  output_fem <- add_line(pheno, "white_euro",  output_fem, pgs_white_euro_fem)
  output_male <- add_line(pheno, "white_euro",  output_male, pgs_white_euro_male)
  #output_all <- add_line(pheno, "afr",  output_all, pgs_afr_all)
  #output_fem <- add_line(pheno, "afr",  output_fem, pgs_afr_fem)
  #output_male <- add_line(pheno, "afr",  output_male, pgs_afr_male)
  #output_all <- add_line(pheno, "asn",  output_all, pgs_asn_all)
  #output_fem <- add_line(pheno, "asn",  output_fem, pgs_asn_fem)
  #output_male <- add_line(pheno, "asn",  output_male, pgs_asn_male)
  
}
head(output_all[, "ci_95_l", ])
save( output_all, output_fem, output_male,file="/Users/manuelacostantino/Documents/descaler_inputs/compiled_output.Rdata")  

