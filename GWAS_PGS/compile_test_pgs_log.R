library("dplyr")

PHENO_DIGITS <- commandArgs(trailingOnly = TRUE)[1]
PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4", ifelse(grepl("FEV1674178", PHENO_DIGITS), "FEV1", ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1", gsub("[0-9]+$", "", PHENO_DIGITS))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)
PHENOS <- "/gpfs/data/ukb-share/extracted_phenotypes/"
TEST_POP <- "whitebrit"
PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")
THRESH <- c("0.0000000001", "0.00000001", "0.000001", "0.0001", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5")
PHENO_FILE <- paste0(PHENOS, PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno")
R2_dirs <- list(
  DEFAULT = "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_default_log_prediction/",
  LOG = "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_log_log_prediction/"
)
PGS_dirs <- list(
  DEFAULT = "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/",
  LOG = "/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out_log/"
)

# Function for R2 calculation
R2fxn <- function(x, y, X) {
  X <- as.matrix(X)
  y <- resid(lm(y ~ X, na.action = na.exclude))
  if (length(x) == length(y)) {
    x <- resid(lm(x ~ X, na.action = na.exclude))
  } else {
    for (j in 1:ncol(x)) {
      xj <- x[, j]
      x[, j] <- resid(lm(xj ~ X, na.action = na.exclude))
    }
  }
  summary(lm(y ~ 1 + x))
}

R2fxn <- function(x, y, X) {
  X <- as.matrix(X)
  y <- resid(lm(y ~ X, na.action = na.exclude))
  if (length(x) == length(y)) {
    x <- resid(lm(x ~ X, na.action = na.exclude))
  } else {
    for (j in 1:ncol(x)) {
      xj <- x[, j]
      x[, j] <- resid(lm(xj ~ X, na.action = na.exclude))
    }
  }
  summary(lm(y ~ 1 + x))
}

pheno_table <- read.table(PHENO_FILE, header = TRUE)
covars <- read.table("/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC674178.pheno", header = TRUE)
covs_table <- covars[, 1:14]
rm(covars)

cat("pheno:", PHENO_LOWER, "\n")

# Loop over scales
for (scale in c("DEFAULT", "LOG")) {

  R2 <- R2_dirs[[scale]]
  PGS_path <- PGS_dirs[[scale]]

  dir.create(R2, showWarnings = FALSE)
  dir.create(paste0(R2, "all_thresh/"), showWarnings = FALSE)

  out_file <- paste0(R2, "processed_", PHENO_LOWER, "_", TEST_POP, "_PGS.txt")
  all_thresh <- paste0(R2, "all_thresh/all_thresh_", PHENO_LOWER, "_", TEST_POP, "_PGS.txt")
  if (file.exists(out_file)) next

  prs.result <- NULL

  for (i in THRESH) try({

    # Subset pheno_cov, y, X for this threshold
    PGS_file <- paste0(PGS_path, PHENO_LOWER, "/", PHENO_LOWER, "_pgs.", i, ".sscore")
    PGS_IDs <- read.table(PGS_file, header = FALSE, col.names = PGS_COLS)[, 1:2]

    covs <- covs_table[covs_table$FID %in% PGS_IDs$FID, ]
    pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ]
    pheno_cov <- merge(pheno, covs, by = c("FID", "IID"))
    colnames(pheno_cov)[3] <- "pheno_code"
    rm(covs, pheno, PGS_IDs)

    # Merge PGS scores
    prs_table <- read.table(PGS_file, header = FALSE, col.names = PGS_COLS)[, c(1,6)]
    prs_df <- merge(pheno_cov, prs_table, by = "FID")
    prs <- prs_df$SCORE1_SUM

    X_input <- prs_df[, !colnames(prs_df) %in% c("FID", "IID", "pheno_code", "SCORE1_SUM")]

    # Compute R2
    if (scale == "DEFAULT") { # different from old script
      prs <- exp(log(prs) )
    }
    prs.R2 <- R2fxn(x = prs, y = log(prs_df$pheno_code), X = X_input)

    prs.result <- rbind(prs.result, data.frame(
      Threshold = i,
      R2 = prs.R2$r.squared,
      w0 = prs.R2$coef["x", "Estimate"],
      type = "pgs"
    ))
  })

  best_std_thresh <- prs.result[which.max(prs.result$R2), ]
  print(best_std_thresh)

  best_prs_table <- read.table(
    paste0(PGS_path, PHENO_LOWER, "/", PHENO_LOWER, "_pgs.", best_std_thresh$Threshold, ".sscore"),
    header = FALSE, col.names = PGS_COLS
  )[, c(1,6)]
  best_prs <- merge(pheno_cov, best_prs_table, by = "FID")$SCORE1_SUM
  rm(best_prs_table)

  write.table(prs.result, file = all_thresh)
  write.table(best_std_thresh, file = out_file)
}