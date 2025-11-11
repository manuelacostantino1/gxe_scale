rm(list=ls())
home_dir <- "/Users/manuelacostantino/Desktop/descaler_analysis/"
load(paste0(home_dir,"compiled_output.Rdata"))

st <- format(Sys.time(), "%m-%d")
pheno_list <- rownames(output_all)
base_cols <- c('#660099', '#EE99AA','#66CCFF','#004488','#DDDDDD')
pop <- "white_euro"

# Clean phenotype names
pheno_clean <- sapply(pheno_list, function(pheno){
  if (grepl("FEV1674178", pheno)) return("FEV1")
  if (grepl("IGF-1674178", pheno)) return("IGF-1")
  if (grepl("EA4", pheno)) return("EA4")
  tolower(gsub("[0-9]+$", "", pheno))
})

plot_dir <- paste0(home_dir,"plots_pgs_r2/")
dir.create(plot_dir, showWarnings = FALSE)

# Set up data
df <- data.frame(
  r2_male = ((output_male[pheno_list,'r20',pop]-output_male[pheno_list,'r21',pop])/output_male[pheno_list,'r21',pop])*100,
  r2_female = ((output_fem[pheno_list,'r20',pop]-output_fem[pheno_list,'r21',pop])/output_fem[pheno_list,'r21',pop])*100,
  CI_25s_male = output_male[pheno_list, 'ci_95_l', pop]*100,
  CI_75s_male = output_male[pheno_list, 'ci_95_u', pop]*100,
  CI_25s_fem = output_fem[pheno_list, 'ci_95_l', pop]*100,
  CI_75s_fem = output_fem[pheno_list, 'ci_95_u', pop]*100,
  pheno_clean = pheno_clean
)

# Sort by female R2
df_sorted <- df[order(df$r2_female), ]
y_pos <- seq_along(df_sorted$pheno_clean)

# Build plot
png(paste0(plot_dir,"r2_sex_",pop,"_vertical.png"),
    width = 10, height = 10, units = 'in', res = 300)
par(mai=c(1,3,1,0.5)) # increase left margin for phenotype labels

# Plot points
plot(df_sorted$r2_male, y_pos, type = "p",
     pch = 16, col = base_cols[3],
     xlab = "R2 % change", ylab = "", yaxt = "n",
     las = 1, cex.lab = 1.5, cex = 1.5)
points(df_sorted$r2_female, y_pos, pch = 16, col = base_cols[1], cex = 1.5)

# Add y-axis phenotype labels
axis(2, at = y_pos, labels = df_sorted$pheno_clean, las = 1, cex.axis = 0.8)

# Vertical zero line
abline(v = 0, col = "black")

# Draw 95% CI as horizontal arrows
ci_cap <- 0.2  # length of end ticks in y-units
segments(x0 = df_sorted$CI_25s_male, x1 = df_sorted$CI_75s_male,
         y0 = y_pos, y1 = y_pos, col = base_cols[3], lwd = 1.5)
segments(x0 = df_sorted$CI_25s_male, x1 = df_sorted$CI_25s_male,
         y0 = y_pos - ci_cap, y1 = y_pos + ci_cap, col = base_cols[3], lwd = 1.5)
segments(x0 = df_sorted$CI_75s_male, x1 = df_sorted$CI_75s_male,
         y0 = y_pos - ci_cap, y1 = y_pos + ci_cap, col = base_cols[3], lwd = 1.5)

segments(x0 = df_sorted$CI_25s_fem, x1 = df_sorted$CI_75s_fem,
         y0 = y_pos, y1 = y_pos, col = base_cols[1], lwd = 1.5)
segments(x0 = df_sorted$CI_25s_fem, x1 = df_sorted$CI_25s_fem,
         y0 = y_pos - ci_cap, y1 = y_pos + ci_cap, col = base_cols[1], lwd = 1.5)
segments(x0 = df_sorted$CI_75s_fem, x1 = df_sorted$CI_75s_fem,
         y0 = y_pos - ci_cap, y1 = y_pos + ci_cap, col = base_cols[1], lwd = 1.5)

# Title and legend
legend("right", legend = c("Males","Females","95% CI"),
       col = c(base_cols[3], base_cols[1], base_cols[3]),
       pch = c(16,16,NA), lty = c(NA,NA,1), lwd=2, cex=1.2, bty="n")

dev.off()

# NOTE THAT I AM INTENTIONALLY FLIPPINF THE SEXES EVE





