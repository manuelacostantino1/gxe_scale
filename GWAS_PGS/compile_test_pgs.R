library("dplyr")

lambda <- 0
PHENO_DIGITS <- commandArgs(trailingOnly = TRUE)[1]
PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4", ifelse(grepl("FEV1674178", PHENO_DIGITS), "FEV1", ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1", gsub("[0-9]+$", "", PHENO_DIGITS))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)
PHENOS <- "/gpfs/data/ukb-share/extracted_phenotypes/"
TEST_POP <- "whitebrit"
PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")
THRESH <- c("0.0000000001", "0.00000001", "0.000001", "0.0001", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5")
PHENO_FILE <- paste0(PHENOS, PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno")
R2_TEST <- "SPEARMAN"

# Crash if R2 is not spearman or pearson
if (!(R2_TEST %in% c("SPEARMAN", "PEARSON"))) {
  stop("R2_TEST must be either 'SPEARMAN' or 'PEARSON'")
}

# Set the correct input and output directories depending on scale and R2 test
if (lambda == 0) {
  if (R2_TEST == "PEARSON") {
    R2 <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_log/"
  } else {
    R2 <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_spearman_log/"
  }
  PGS <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out_log/"
} else if (lambda == 1) {
  if (R2_TEST == "PEARSON") {
    R2 <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out/"
  } else {
    R2 <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_spearman/"
  }
  PGS <- "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/"
}


R2fxn <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid( lm( y ~ X ,na.action = na.exclude) )
  if( length(x) == length(y) ){
    x <- resid( lm( x ~ X ,na.action = na.exclude) )
  } else {
    for( j in 1:ncol(x) ){
      xj <- x[,j]
      x[,j] <- resid( lm( xj ~ X ,na.action = na.exclude) )
    }
  }
  summary( lm( y ~ 1 + x ) )
}

R2fxn_spearman <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid( lm( y ~ X ,na.action = na.exclude) )
  if( length(x) == length(y) ){
    x <- resid( lm( x ~ X ,na.action = na.exclude) )
  } else {
    for( j in 1:ncol(x) ){
      xj <- x[,j]
      x[,j] <- resid( lm( xj ~ X ,na.action = na.exclude) )
    }
  }
  (cor(x, y, method = "spearman", use = "complete.obs"))^2
}

pheno_table <- read.table(PHENO_FILE, header=TRUE)

cat("pheno:", PHENO_LOWER, "\n")

dir.create(R2, showWarnings = FALSE)
dir.create(paste0(R2,"all_thresh/"), showWarnings = FALSE)
out_file <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_PGS.txt")
all_thresh <- paste0(R2,"all_thresh/all_thresh_",PHENO_LOWER,"_",TEST_POP,"_PGS.txt")
if(file.exists(out_file)) {
  quit(save = "no", status = 0)
}

### build phenos+covariates for PGS individuals
covars <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC674178.pheno",header=TRUE) 
covs_table <- covars[, 1:14]
rm(covars)
PGS_file <- paste0(PGS,PHENO_LOWER,"/", PHENO_LOWER,"_pgs.",THRESH[1],".sscore")  # This file only used to grab IDs
PGS_IDs <- read.table(PGS_file,header=FALSE, col.names=PGS_COLS)[,1:2]
covs <- covs_table[covs_table$FID %in% PGS_IDs$FID, ] # only keep the people in the prs
pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ] # only keep the people in the prs
pheno_cov <- merge(pheno, covs, by=c("FID", "IID"))
colnames(pheno_cov)[3] <- "pheno_code"
rm( PGS_IDs,PGS_file,covs_table,pheno )
y <- pheno_cov$pheno_code
X <- pheno_cov[,!colnames(pheno_cov)%in%c("FID","IID","pheno_code")] 
  
### Get the results for all threshold
prs.result <- NULL	
for(i in THRESH){

  prs_table <- read.table(paste0(PGS, PHENO_LOWER, "/", PHENO_LOWER, "_pgs.", i, ".sscore"),
                          header = FALSE, col.names = PGS_COLS)[, c(1,6)]
  prs <- merge(pheno_cov, prs_table, by = "FID")$SCORE1_SUM
  
  if (lambda == 0) {
    prs <- exp(prs)
  }

  if (R2_TEST == "SPEARMAN") {
    prs.R2 <- R2fxn_spearman(x = prs, y = y, X = X)
    print(prs.R2)
    prs.result <- rbind(prs.result, data.frame(Threshold = i, R2 = prs.R2, w0 = NA, type = "pgs"))
  } else {
    mod <- R2fxn(x = prs, y = y, X = X)
    prs.result <- rbind(prs.result, data.frame(Threshold = i, R2 = mod$r.squared, w0 = mod$coef["x", "Estimate"], type = "pgs"))
    print(prs.result)
  }
}
best_std_thresh <- prs.result[which.max(prs.result$R2),]
print(best_std_thresh)
rm(prs_table)

best_prs_table <- read.table(paste0(PGS,PHENO_LOWER,"/", PHENO_LOWER,"_pgs.",best_std_thresh[[1]],".sscore"), header=FALSE, col.names=PGS_COLS)[,c(1,6)]
best_prs <- merge(pheno_cov,best_prs_table,by="FID")$SCORE1_SUM
rm(best_prs_table)

prs_full <- prs.result
prs.update <- best_std_thresh
  
write.table(prs_full, file = all_thresh)
write.table(prs.update, file = out_file)
