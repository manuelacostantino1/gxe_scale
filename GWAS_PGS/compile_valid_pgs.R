source("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/descaler/scripts/bootstrap_prs.R")
library("dplyr") 

######
R2_TEST <-"SPEARMAN"
PHENO_DIGITS <- commandArgs(trailingOnly = TRUE)[1]
PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4", ifelse(grepl("FEV1674206", PHENO_DIGITS), "FEV1", ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1", gsub("[0-9]+$", "", PHENO_DIGITS))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)
PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")
PHENO_FILE <- paste0("/gpfs/data/ukb-share/extracted_phenotypes/", PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno") 
pheno_table <- read.table(PHENO_FILE, header=TRUE)
######

if (!(R2_TEST %in% c("SPEARMAN", "PEARSON"))) {
  stop("R2_TEST must be either 'SPEARMAN' or 'PEARSON'")
}

if (R2_TEST == "SPEARMAN") {
  R20   <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_spearman_log/"
} else{
  R20   <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_log/"
}
best_threshes <- read.table(paste0(R20, "processed_", PHENO_LOWER, "_whitebrit_PGS.txt"), header = TRUE)
thresh0   <- format(best_threshes[1,1],scientific=F) 

if (R2_TEST == "SPEARMAN") {
  R21   <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_spearman/"
} else{
  R21   <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out/"
}
best_threshes <- read.table(paste0(R21, "processed_", PHENO_LOWER, "_whitebrit_PGS.txt"), header = TRUE)
thresh1 <- format(best_threshes[1,1],scientific=F)

cat("pheno:", PHENO_LOWER, "\n")
for( group in c( 'all', 'female', 'male' ) ){
for(valid_pop in c("white_euro","afr","asn")){
  if (R2_TEST == "SPEARMAN") {
    out_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs_spearman/log_vs_lin_",PHENO_LOWER,"_",valid_pop,'_', group, ".txt")
    PGS0_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out_log/valid_spearman/",PHENO_LOWER,"_",valid_pop,"_pgs_valid.",thresh0,".sscore")
    PGS1_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/valid_spearman/",PHENO_LOWER,"_",valid_pop,"_pgs_valid.",thresh1,".sscore")
  } else{
    out_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs_noexp/log_vs_lin_",PHENO_LOWER,"_",valid_pop,'_', group, ".txt")
    PGS0_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out_log/valid/",PHENO_LOWER,"_",valid_pop,"_pgs_valid.",thresh0,".sscore")
    PGS1_file <- paste0("/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/valid/",PHENO_LOWER,"_",valid_pop,"_pgs_valid.",thresh1,".sscore")
  }
  print(out_file)
  if(file.exists(out_file)) next()
  
  ### build phenos+covariates and subset to PGS individuals
  covs_table <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC674178.pheno",header=TRUE)[,1:14]
  PGS_IDs <- read.table(PGS0_file,header=FALSE, col.names=PGS_COLS)[,1:2]
  covs <- covs_table[covs_table$FID %in% PGS_IDs$FID, ] # only keep the people in the prs
  if( group == 'female' ){
    covs <- covs[covs['X31.0.0']==1, ]
  } else if( group == 'male' ){
    covs <- covs[covs['X31.0.0']==0, ]
  }

  pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ] # only keep the people in the prs
  pheno_cov <- merge(pheno, covs, by=c("FID", "IID"))
  colnames(pheno_cov)[3] <- "pheno_code"
  rm( PGS_IDs,covs,covs_table,pheno )
  
  y <- pheno_cov$pheno_code
  X <- pheno_cov[,!colnames(pheno_cov)%in%c("FID","IID","pheno_code")]

  ### Step 1: Score all of the tau0
  prs0 <- read.table(PGS0_file, header=FALSE,col.names=PGS_COLS)[,c(1,6)]
  prs1 <- read.table(PGS1_file, header=FALSE,col.names=PGS_COLS)[,c(1,6)]
  prs0 <- exp(merge(pheno_cov,prs0,by="FID")$SCORE1_SUM)
  prs1 <- merge(pheno_cov,prs1,by="FID")$SCORE1_SUM
  if (R2_TEST == "SPEARMAN") {
    prs.R2 <- r2_diff_spearman(prs0,prs1,y=y,X=X,B=1000)
  } else{
    prs.R2 <- r2_diff_boot(prs0,prs1,y=y,X=X,B=1000)
  }
  prs.result <- data.frame(thresh0=thresh0, thresh1=thresh1, r20=prs.R2$r20, r21=prs.R2$r21, ci_95=prs.R2$ci_95, pval=prs.R2$pval )

  write.table(prs.result, file = out_file)

}
}

