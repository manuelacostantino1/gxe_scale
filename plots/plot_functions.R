prs_list <- c("BMI", "age_menopause","chrons", "melanoma", "alzheimers", "coeliac", "osteoporosis", "ulcerative_colitis", "HBA1C", "arthritis", "coronary_artery", "ovarian_cancer", "venous_thromboembolic", "HDL", "asthma", "glaucoma", "parkinsons", "LDL", "atrial_fibrillation", "height", "prostate_cancer", "MS", "bipolar", "hypertension", "psoriasis",  "bone_mineral_density", "intraocular_pressure", "schizophrenia", "T1D", "bowel_cancer", "macular_degenration", "lupus", "stroke", "T2D", "breast_cancer","CVD")
pheno_list <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
  "Calcium", "SystolicBP_auto", "Height", "Testosterone",  "Creatinine", "Triglycerides", 
   "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
  "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol",  "Hip_circumference",
  "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin")
env_list <- c("sex", "age", "Alcohol_intake_frequency", "Statins", "BMI", "EA4", "Smoking_status")
lambdaseq <- c(-3.00000000, -2.87755102, -2.75510204, -2.63265306, -2.51020408, -2.38775510, -2.26530612, -2.14285714, -2.02040816, -1.89795918, -1.77551020, -1.65306122, -1.53061224, -1.40816327, -1.28571429, -1.16326531, -1.04081633, -0.91836735, -0.79591837, -0.67346939, -0.55102041, -0.42857143, -0.30612245, -0.18367347, -0.06122449,  0.00000000,  0.06122449,  0.18367347, 0.30612245,  0.42857143,  0.55102041,  0.67346939,  0.79591837,  0.91836735,  1.00000000, 1.04081633,  1.16326531,  1.28571429,  1.40816327,  1.53061224,  1.65306122,  1.77551020, 1.89795918,  2.02040816)
#lambdashort <- c(0.00000000,  0.06122449,  0.18367347, 0.30612245,  0.42857143,  0.55102041,  0.67346939,  0.79591837,  0.91836735,  1.00000000, 1.04081633,  1.16326531,  1.28571429,  1.40816327,  1.53061224,  1.65306122,  1.77551020, 1.89795918,  2.02040816)
lambdashort <- c(-1.04081633, -0.91836735, -0.79591837, -0.67346939, -0.55102041, -0.42857143, -0.30612245, -0.18367347, -0.06122449,  0.00000000,  0.06122449,  0.18367347, 0.30612245,  0.42857143,  0.55102041,  0.67346939,  0.79591837,  0.91836735,  1.00000000, 1.04081633,  1.16326531,  1.28571429,  1.40816327,  1.53061224,  1.65306122,  1.77551020, 1.89795918,  2.02040816)
custom_palette <- c("Sex" = "#00CCFF", "Age" = "#D55E00", "Alcohol intake frequency" = "#009E73", "Statins" = "#990033", "BMI" = "#FF99FF", "EA4" = "#FFCC00", "Smoking status" = "#003399", "Geometric mean" = "#000000")
threshold <- 0.05/(length(prs_list))
#threshold <- 0.05

make_3_arrays <- function(big_array, threshold, lambdaseq){
  print("slay")
  cat("slay\n")
  # Subset to correct indices
  lambda_indices <- match(dimnames(big_array)[[4]], lambdaseq)
  lambda_indices <- which(!is.na(lambda_indices))
  big_array <- big_array[,,,lambda_indices]
  dimnames(big_array)[[4]] <- as.character(lambdaseq)
  
  unkillable_array <- big_array
  killable_array <- big_array
  null_array <- big_array
  for (i in 1:length(pheno_list)) {
    for (j in 1:length(prs_list)) {
      for (k in 1:length(env_list)) {
        if (all(is.na(big_array[j, i, k, ]))){
          next()
        }
        if (all(big_array[j, i, k, ] < threshold, na.rm=TRUE)) {
          null_array[j, i, k, ] <- NA
          killable_array[j, i, k, ] <- NA
        } else if (big_array[j, i, k, as.character(lambdaseq) == "1"] < threshold){
          unkillable_array[j, i, k, ] <- NA
          null_array[j, i, k, ] <- NA
        } else {
          unkillable_array[j, i, k, ] <- NA
          killable_array[j, i, k, ] <- NA
        }
      }
    }
  }
  return(list(unkillable = unkillable_array,killable = killable_array,null = null_array))
}


make_3_arrays_v2 <- function(big_array, threshold, lambdaseq){

  # Subset to correct indices
  lambda_indices <- match(dimnames(big_array)[[4]], lambdaseq)
  lambda_indices <- which(!is.na(lambda_indices))
  big_array <- big_array[,,,lambda_indices]
  dimnames(big_array)[[4]] <- as.character(lambdaseq)
  
  unkillable_array <- big_array
  killable_array <- big_array
  null_array <- big_array
  for (i in 1:length(pheno_list)) {
    for (j in 1:length(prs_list)) {
      for (k in 1:length(env_list)) {
        if (all(is.na(big_array[j, i, k, ]))){
          next()
        }
        if ((big_array[j, i, k, as.character(lambdaseq) == "1"] < threshold) & (big_array[j, i, k, as.character(lambdaseq) == "0"] < threshold)) {
          null_array[j, i, k, ] <- NA
          killable_array[j, i, k, ] <- NA
        } else if (big_array[j, i, k, as.character(lambdaseq) == "1"] < threshold){
          unkillable_array[j, i, k, ] <- NA
          null_array[j, i, k, ] <- NA
        } else {
          unkillable_array[j, i, k, ] <- NA
          killable_array[j, i, k, ] <- NA
        }
      }
    }
  }
  return(list(unkillable = unkillable_array,killable = killable_array,null = null_array))
}


collapse_interactions <- function(input_array, envs){
  input_array <- input_array[,,,1] 
  df <- data.frame(matrix(0, ncol = length(pheno_list), nrow = (length(envs) * length(prs_list))))
  colnames(df) <- pheno_list  
  i <- 1
  for (env in envs){
    for (prs in prs_list){
      
      interaction <- paste(prs, "_PRS_by_", env, sep="")  
      non_na_indices <- which(!is.na(input_array[prs, , env]))
      df[i, non_na_indices] <- 1  
      rownames(df)[i] <- interaction  
      i <- i + 1  
    }
  }
  return(df)  
}

capitalize_first <- function(s) {
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

lower_first <- function(s) {
  paste0(tolower(substring(s, 1, 1)), substring(s, 2))
}

format_pgs_label <- function(s) {
  s <- gsub("_", " ", s)
  s <- capitalize_first(s)
  paste0(s, " PGS")
}

make_long <- function(input_df){
  input_df$interaction <- rownames(input_df)
  input_df_long <- pivot_longer(input_df, cols = -interaction, names_to = "phenotype", values_to = "Category")
  input_df_long$Category <- as.factor(input_df_long$Category)
  input_df_long$interaction <- gsub("_PRS_by_.*", "", input_df_long$interaction)
  input_df_long$interaction <- gsub("_", " ", input_df_long$interaction)
  input_df_long$interaction <- capitalize_first(input_df_long$interaction)
  input_df_long$phenotype <- gsub("_", " ", input_df_long$phenotype)
  return(input_df_long)
}



