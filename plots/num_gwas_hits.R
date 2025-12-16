library(data.table)
library(GenomicRanges)
library(ggplot2)
library(VennDiagram)
library(grid)

all_phenos <- c("Birth_weight", "HDL", "Alcohol_intake_frequency", "C-reactive_protein", "HbA1c", "Mean_corpuscular_volume",
    "Calcium", "SystolicBP_auto", "Height", "Testosterone", "Arm_fat-free_mass_left", "Creatinine", "Triglycerides", 
    "Arm_fat-free_mass_right", "DiastolicBP_auto", "Platelet_count", "Urate", "Pulse_rate",
    "Urea", "FEV1", "IGF-1", "RBC", "Vitamin_D",  "LDL", "Leukocyte_count", "SHBG", "BMI", "EA4", "Cholesterol", "Hip_circumference",
    "Arm_fat-free_mass_avg", "Hip_to_waist", "Waist_circumference", "Glucose", "Whole_body_fat_mass", "WHRadjBMI_Zhu")

df <- data.frame(
  Phenotype = character(),
  Shared_default = integer(),
  Specific_default = integer(),
  Shared_log = integer(),
  Specific_log = integer()
)

for (pheno in all_phenos) {
  pheno_lower <- tolower(pheno)
  cat("Processing:", pheno_lower, "\n")
  
  tryCatch({
    merge_dist <- 500000  
    comp_dist  <- 250000 

    # Load files
    clumps_default <- read.table(
      paste0("/scratch/mcostantino/pgs_output/clumped_genos/ukb_chr1-22_whitebrit_", pheno_lower, ".clumps"),
      header = TRUE, comment.char = "", check.names = FALSE
    )
    clumps_log <- read.table(
      paste0("/scratch/mcostantino/pgs_output/clumped_genos_log/ukb_chr1-22_whitebrit_", pheno_lower, ".clumps"),
      header = TRUE, comment.char = "", check.names = FALSE
    )

    # Filter by significance
    clumps_default <- clumps_default[clumps_default$P < 5e-8, ]
    clumps_log     <- clumps_log[clumps_log$P < 5e-8, ]

    # Merge hits
    merge_hits <- function(df, merge_dist) {
      gr <- GRanges(
        seqnames = paste0("chr", df$`#CHROM`),
        ranges = IRanges(start = df$POS, end = df$POS),
        SNP = df$ID, P = df$P
      )
      merged <- reduce(gr, min.gapwidth = merge_dist + 1)
      keep_snps <- sapply(seq_along(merged), function(i) {
        overlaps <- findOverlaps(gr, merged[i])
        region_snps <- gr[queryHits(overlaps)]
        region_snps[which.min(region_snps$P)]$SNP
      })
      df[df$ID %in% keep_snps, ]
    }

    clumps_default_merged <- merge_hits(clumps_default, merge_dist)
    clumps_log_merged     <- merge_hits(clumps_log, merge_dist)

    # Create GRanges for ±500kb
    make_loci <- function(df) {
      GRanges(
        seqnames = paste0("chr", df$`#CHROM`),
        ranges = IRanges(
          start = pmax(df$POS - comp_dist, 1),
          end   = df$POS + comp_dist
        ),
        SNP = df$ID
      )
    }

    loci1 <- make_loci(clumps_default_merged)
    loci2 <- make_loci(clumps_log_merged)

    overlap <- findOverlaps(loci1, loci2)
    shared_default <- length(unique(queryHits(overlap)))
    shared_log     <- length(unique(subjectHits(overlap)))
    specific_default <- length(loci1) - shared_default
    specific_log     <- length(loci2) - shared_log

    df <- rbind(df, data.frame(
      Phenotype = pheno,
      Shared_default = shared_default,
      Specific_default = specific_default,
      Shared_log = shared_log,
      Specific_log = specific_log
    ))
    print(df)
  }, error = function(e) {
      cat("Error processing", pheno_lower, ":", e$message, "\n")
       })
}

df