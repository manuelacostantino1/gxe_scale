# Load required packages
library(dplyr)

# Set phenotype and paths
pheno <- "testosterone"
pheno_title <- gsub("_", " ", pheno)

lambda <- 1
print(paste("lambda=", lambda))

GWAS <- "/gpfs/data/ukb-share/dahl/manuela/gwas_results/"
TRAIN_POP <- "whitebrit"
GWAS_COLS <- c("CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", 
               "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", 
               "P", "ERRCODE")

# Update paths based on lambda
if (lambda == 0) {
  JOINED_GWAS <- "/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/"
} else if (lambda == 1) {
  JOINED_GWAS <- paste("/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/", sep = "")
}

# Load GWAS results
gwas_file <- paste0(JOINED_GWAS, pheno, "/chr1-22_", TRAIN_POP, "_", pheno, ".assoc.linear")
gwas_data <- read.table(gwas_file, col.names = GWAS_COLS)

# Calculate chi-square statistics from p-values
gwas_data <- gwas_data %>% 
  filter(!is.na(P) & P > 0)  # Exclude NA or zero p-values
gwas_data <- gwas_data %>% 
  mutate(chi2_stat = qchisq(1 - P, df = 1))

# Calculate genomic control factor (lambda GC)
lambda_gc <- median(gwas_data$chi2_stat, na.rm = TRUE) / qchisq(0.5, df = 1)
print(paste("Genomic Control Factor (lambda GC):", lambda_gc))

# Save results
gc_result_path <- paste0(JOINED_GWAS, pheno, "_lambda_gc.txt")
write(paste("Phenotype:", pheno, "Lambda GC:", lambda_gc), gc_result_path)

print("Lambda GC calculation complete and saved.")