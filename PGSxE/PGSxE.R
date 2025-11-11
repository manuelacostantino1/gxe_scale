# Take in command line arguments 
taskid <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
env <- commandArgs(trailingOnly = TRUE)[2]
my_pheno <- commandArgs(trailingOnly = TRUE)[3]
my_prs <- commandArgs(trailingOnly = TRUE)[4]

outliers <- TRUE
error_log_file <- "/scratch/mcostantino/descaler/logs/error_log.txt"

# Make list of interesting phenotypes and prs
all_phenos <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Arm_fat-free_mass_left", "Creatinine", "Triglycerides", 
    "Arm_fat-free_mass_right", "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D", "Basophill_count", "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol", "Eosinophill_count", "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin")
all_prs <- c("sPRS_BMI", "sPRS_age_menopause", "sPRS_chrons", "sPRS_melanoma", "sPRS_CVD", "sPRS_alzheimers", "sPRS_coeliac", "sPRS_osteoporosis", "sPRS_ulcerative_colitis",
	     "sPRS_HBA1C", "sPRS_arthritis", "sPRS_coronary_artery", "sPRS_ovarian_cancer", "sPRS_venous_thromboembolic", "sPRS_HDL", "sPRS_asthma", "sPRS_glaucoma", "sPRS_parkinsons", "sPRS_LDL", 
	     "sPRS_atrial_fibrillation", "sPRS_height", "sPRS_prostate_cancer", "sPRS_MS", "sPRS_bipolar", "sPRS_hypertension", "sPRS_psoriasis", "sPRS_bone_mineral_density", "sPRS_intraocular_pressure",
	     "sPRS_schizophrenia", "sPRS_T1D", "sPRS_bowel_cancer", "sPRS_lupus", "sPRS_stroke", "sPRS_T2D", "sPRS_breast_cancer", "sPRS_macular_degenration")
tryCatch({
# Figure out what the phenos you want are
x <- as.matrix(expand.grid( 1:length(all_phenos), 1:length(all_prs) )[taskid,])
my_pheno <- all_phenos[x[1]]
my_prs <- all_prs[x[2]]
print(paste("Phenotype being used: ", my_pheno))
print(paste("PRS being used: ", my_prs))
write(paste("Taskid: ", taskid, " Phenotype: ", my_pheno, " PRS: ", my_prs, " env: ", env), file = "/scratch/mcostantino/descaler/logs/id_pheno_prs.txt", append = TRUE)

# Check environment
if (env %in% c("sex","age","EA4","BMI", "Smoking_status", "Alcohol_intake_frequency", "Statins")){
    print(paste("Environment being used: ", env))
} else {
    print(paste("Environment", env, "not recognized"))
    stop()
}

# Path to ouptut directory
if (outliers ==FALSE) {
  output_dir <- paste("/gpfs/data/ukb-share/dahl/manuela/scale_project/ukb_prs_results_", env, sep="")
} else {
  output_dir <- paste("/gpfs/data/ukb-share/dahl/manuela/PGSxE/ukb_prs_results_", env, sep="")
}
#output_dir <- paste("/scratch/mcostantino/descaler/test_gxe/ukb_prs_results_", env, sep="")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
subdirs <- c("pvals", "stder", "r2", "gmain", "emain", "coefs")
for (subdir in subdirs) {
  if (env != "Smoking_status" || !(subdir %in% c("gmain", "emain", "coefs"))) {
    dir.create(file.path(output_dir, paste0(subdir, "_", env)), recursive = TRUE, showWarnings = FALSE)
  }
}


expected_files <- c(
  paste0(output_dir, "/pvals_", env, "/pvals_", my_pheno, "_", my_prs, ".txt"),
  paste0(output_dir, "/stder_", env, "/stder_", my_pheno, "_", my_prs, ".txt"),
  paste0(output_dir, "/r2_", env, "/r2_", my_pheno, "_", my_prs, ".txt")
)

# For environments other than Smoking_status, add the other expected files
if (env != "Smoking_status") {
  expected_files <- c(expected_files,
    paste0(output_dir, "/gmain_", env, "/gmain_", my_pheno, "_", my_prs, ".txt"),
    paste0(output_dir, "/emain_", env, "/emain_", my_pheno, "_", my_prs, ".txt"),
    paste0(output_dir, "/coefs_", env, "/coefs_", my_pheno, "_", my_prs, ".txt")
  )
}

# Skip this job if all expected output files already exist
if (all(file.exists(expected_files))) {
  message("All output files for this pheno-PRS pair already exist. Skipping job.")
  quit(save = "no", status = 0)
}


# Open the correct files
prs_path <- paste("/gpfs/data/ukb-share/dahl/manuela/ukb_prs_downloading/polygenic_scores/", my_prs, sep="")
prs_file <- list.files(prs_path, full.names = TRUE)[1]
prs_table <- read.table(prs_file, header = TRUE)
pheno_path <- paste("/gpfs/data/ukb-share/extracted_phenotypes/", my_pheno, sep="")
pheno_file <- list.files(pheno_path, full.names = TRUE)[1]
pheno_table <- read.table(pheno_file, header = TRUE)
cov_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covar_full/covar_full_age2.pheno", header=TRUE)
e_path <- paste("/gpfs/data/ukb-share/extracted_phenotypes/",env,sep="")
merged_data <- merge(pheno_table, prs_table, by = "IID")
merged_data <- merge(merged_data, cov_table, by = "IID")
if (!(env %in% c("age", "sex"))){
    e_file <- list.files(path = e_path, full.names = TRUE)[1]
    e_table <- read.table(e_file, header=TRUE)
    if (env == "Statins"){
        e_table <- e_table[,1:3]
    }
    merged_data <- merge(merged_data, e_table, by = "IID")
}

# Subset to whitebrits
famfile <- read.table("/scratch/reneefonseca/genotypes/whitebrit/ukb_chr10_whitebrit_QC.fam")
colnames(famfile) <- c("FID", "IID", "PID", "MID", "sex", "phen")
merged_data <- merge(merged_data, famfile, by = "IID")

# Build data frame with prs, pheno, covariates (except the pcs)
if (env %in% c("Smoking_status","Statins")) {
  df <- data.frame(pheno = merged_data[,3], prs = merged_data[,5], env = as.factor(merged_data[,(ncol(merged_data)-5)]), age = merged_data$X21003.0.0, age2 = (merged_data$X21003.0.0)^2, center = as.factor(merged_data$X54.0.0), sex = as.factor(merged_data$X31.0.0))
} else if (env == "age"){
  df <- data.frame(pheno = merged_data[,3], prs = merged_data[,5], env = merged_data$X21003.0.0, center = as.factor(merged_data$X54.0.0), sex = as.factor(merged_data$X31.0.0))
} else if (env == "sex") {
  df <- data.frame(pheno = merged_data[,3], prs = merged_data[,5], env = as.factor(merged_data$X31.0.0), center = as.factor(merged_data$X54.0.0), age = merged_data$X21003.0.0, age2 = merged_data$age2)
} else {
  df <- data.frame(pheno = merged_data[,3], prs = merged_data[,5], env = merged_data[,(ncol(merged_data)-5)], age = merged_data$X21003.0.0, age2 = (merged_data$X21003.0.0)^2, center = as.factor(merged_data$X54.0.0), sex = as.factor(merged_data$X31.0.0))
}

# Add PCs
pc1_index <- which(colnames(merged_data)=="X22009.0.1")
pcs <- merged_data[pc1_index:(pc1_index+39)]
colnames(pcs) <- paste("PC", 1:40, sep = "")
df <- cbind(df,pcs)

# Remove people with missing data and scale, remove outliers if needed
quants <- quantile(df$pheno, c(0.001, 0.999), na.rm=TRUE)
if (outliers == TRUE){
  df$pheno[df$pheno < quants[1]] <- NA
  df$pheno[df$pheno > quants[2]] <- NA
}
df <- na.omit(df)
df$prs <- scale(df$prs)

### Probs not needed anymore? from when there were negative values in some phenos
print(env)
if (env == "Smoking_status"){
  df <- df[as.numeric(as.character(df$env))>=0,]
  df$env <- droplevels(df$env)
}
df <- df[df$pheno>0,]  # this for sure keep

# Define boxcox function 
boxcox <- function(y,lambda){
  if( lambda == 0 ) return(scale(log(y)))
  y <-(y^lambda-1) / lambda
  y <- scale(y)
}

# Run a linear model and get the regression coefficient and p values of the interaction
lambdas <- c(-3.00000000, -2.87755102, -2.75510204, -2.63265306, -2.51020408, -2.38775510, -2.26530612, -2.14285714, -2.02040816, -1.89795918, -1.77551020, -1.65306122, -1.53061224, -1.40816327, -1.28571429, -1.16326531, -1.04081633, -0.91836735, -0.79591837, -0.67346939, -0.55102041, -0.42857143, -0.30612245, -0.18367347, -0.06122449,  0.00000000,  0.06122449,  0.18367347, 0.30612245,  0.42857143,  0.55102041,  0.67346939,  0.79591837,  0.91836735,  1.00000000, 1.04081633,  1.16326531,  1.28571429,  1.40816327,  1.53061224,  1.65306122,  1.77551020, 1.89795918,  2.02040816,  2.14285714,  2.26530612,  2.38775510,  2.51020408,  2.63265306,2.75510204,  2.87755102,  3.00000000)
#lambdas <- c(0, 1)

descaler <- function( g, e, y, covars, lambdaseq=lambdas ){
  coefs <- c()
  ps <- c()
  stder <- c()
  r2 <- c()
  gmain <- c()
  emain <- c()
  covar_names <- colnames(covars)
  for( j in seq(along=lambdaseq) ){
    print(j)
    y1 <- boxcox(y,lambdaseq[j])
    print(head(y1))
    interaction_terms <- paste("e*covars$", covar_names, sep = "")
    formula <- as.formula(paste("y1 ~ g * e +", paste("covars$",covar_names, collapse = " + "), "+", paste(interaction_terms, collapse = " + ")))
    mod <- lm(formula)
    r2[j] <- summary(mod)$r.squared
    if (env == "Smoking_status"){
        print("treating as multi-level factor")
        formula <- as.formula(paste("y1 ~ g + e +", paste("covars$",covar_names, collapse = " + "), "+", paste(interaction_terms, collapse = " + ")))
        mod2 <- lm(formula)
        anova_result <- anova(mod, mod2)
        ps[j] <- anova_result[2, "Pr(>F)"]
        stder[j] <- summary(mod)$coef[2,"Std. Error"]
        # ADD RETURN MAIN EFFECTS and COEFFICIENTS
    } else if (env %in% c("sex", "Statins")){
      print("Treating as binary factor")
          coefs[j] <- summary(mod)$coef['g:e1',"Estimate"]
          ps[j] <- summary(mod)$coef['g:e1',"Pr(>|t|)"]
          stder[j] <- summary(mod)$coef['g:e1',"Std. Error"]
          gmain[j] <- summary(mod)$coef['g',"Pr(>|t|)"]
          emain[j] <- summary(mod)$coef['e1',"Pr(>|t|)"]
    } else {
      print("not treating as factor")
          coefs[j] <- summary(mod)$coef['g:e',"Estimate"]
          ps[j] <- summary(mod)$coef['g:e',"Pr(>|t|)"]
          stder[j] <- summary(mod)$coef['g:e',"Std. Error"]
          gmain[j] <- summary(mod)$coef['g',"Pr(>|t|)"]
          emain[j] <- summary(mod)$coef['e',"Pr(>|t|)"]
    }
  }
  return( list( coefs=coefs, pvals=ps, stder=stder, r2=r2 , gmain=gmain, emain=emain) ) 
}

results <- descaler(df$prs, df$env, df$pheno, df[4:ncol(df)])
write.table(results$pvals, file.path(paste(output_dir,"/pvals_", env,"/pvals_",my_pheno,"_",my_prs,".txt", sep="")))
write.table(results$stder, file.path(paste(output_dir,"/stder_", env,"/stder_",my_pheno,"_",my_prs,".txt", sep="")))
write.table(results$r2, file.path(paste(output_dir,"/r2_", env,"/r2_",my_pheno,"_",my_prs,".txt", sep="")))
if (env != "Smoking_status"){
    write.table(results$gmain, file.path(paste(output_dir,"/gmain_", env,"/gmain_",my_pheno,"_",my_prs,".txt", sep="")))
    write.table(results$emain, file.path(paste(output_dir,"/emain_", env,"/emain_",my_pheno,"_",my_prs,".txt", sep="")))
    write.table(results$coefs, file.path(paste(output_dir,"/coefs_", env,"/coefs_",my_pheno,"_",my_prs,".txt", sep="")))
}

}, error = function(e) {
    # Log the error details to the error log file
    write(paste("Taskid: ", taskid, " Phenotype: ", my_pheno, " PRS: ", my_prs, " env: ", env, 
                " Error: ", e$message), file = error_log_file, append = TRUE)
    print(paste("Error encountered. Logged to:", error_log_file))
})