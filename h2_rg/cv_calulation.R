# Make list of phenotypes and make vectors
phenotypes <- c("Alcohol_intake_frequency674178", "Arm_fat-free_mass_left674178", "Basophill_count674178", "Birth_weight674178", "BMI674178", "Calcium674178", "Cholesterol674178", "Creatinine674178", "DiastolicBP_auto674178", "Eosinophill_count674178", "Glucose674178", "HDL674178", "HbA1c674178", "Height674178", "Hip_circumference674178", "IGF-1674178", "LDL674178", "Leukocyte_count674178", "Mean_corpuscular_volume674178", "Platelet_count674178", "Pulse_rate674178", "RBC674178", "SHBG674178", "SystolicBP_auto674178", "Testosterone674178", "Triglycerides674178", "Urate674178", "Urea674178", "Vitamin_D674178", "Waist_circumference674178", "Whole_body_fat_mass674178", "WHRadjBMI_Zhu", "WHRadjBMI_Emdin", "FEV1674178", "EA4")
cv <- c()
delta.cv <- c()
phenos <- c()

# Open file with sex
covars <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC674178.pheno",header=TRUE) 

# Iterate through phenotypes
for (pheno in phenotypes){

    # Open corresponding file
    pheno_no_digits <- ifelse(grepl("EA4", pheno), "EA4", ifelse(grepl("FEV1674178", pheno), "FEV1", ifelse(grepl("IGF-1674178", pheno), "IGF-1", gsub("[0-9]+$", "", pheno))))
    
    path <- paste("/gpfs/data/ukb-share/extracted_phenotypes/", pheno_no_digits, "/", pheno, ".pheno", sep="")
    pheno_data <- read.table(path, header=TRUE)
    print(pheno)
    # Calculate CV
    coefficient <- (sd(pheno_data[,3], na.rm=TRUE)/mean(pheno_data[,3],na.rm=TRUE))
    cv <- c(cv, coefficient)

    # Calculate delta CV
    big_df <- merge(pheno_data, covars, by=c("FID", "IID"))
    male_pheno <- big_df[which(big_df$X31.0.0==1),3]
    female_pheno <- big_df[which(big_df$X31.0.0==0),3]
    cv_male <- sd(male_pheno, na.rm=TRUE)/mean(male_pheno, na.rm=TRUE)
    print(cv_male)
    print(sd(male_pheno, na.rm=TRUE))
    print(mean(male_pheno, na.rm=TRUE))
    cv_female <- sd(female_pheno, na.rm=TRUE)/mean(female_pheno, na.rm=TRUE)
    delta.cv <- c(delta.cv, (cv_male-cv_female))
    print(delta.cv[length(delta.cv)])

    phenos <- c(phenos, tolower(pheno_no_digits))
}

df <- data.frame(phenos=phenos, cv=cv, dcv=delta.cv)
write.table(df,"/scratch/mcostantino/descaler/coefficient_variation.txt")