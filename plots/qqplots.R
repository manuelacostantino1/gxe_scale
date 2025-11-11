# https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/
# https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html

# load packages
library(dplyr)
library(ggplot2)  

pheno="height"
pheno_title <- gsub("_", " ", pheno)

lambda=0
lambda_no_dots=0
print(paste("lambda=",lambda))

GENOS="/scratch/reneefonseca/genotypes/"
PHENOS="/gpfs/data/ukb-share/extracted_phenotypes/"
GWAS="/gpfs/data/ukb-share/dahl/manuela/gwas_results/"
LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
TRAIN_POP="whitebrit"
JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/"
GWAS_COLS <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF?","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
QQ_PLOTS="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/gwas/"

if(lambda==0){
  print(paste("Lambda=", lambda))
  GWAS="/gpfs/data/ukb-share/dahl/manuela/gwas_results_log/"
	JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/"
  QQ_PLOTS="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/gwas_log/"
} else if (lambda!= 1) {
  print(paste("Lambda=", lambda))
  GWAS=paste("/gpfs/data/ukb-share/dahl/manuela/gwas_results_",lambda_no_dots,"/",sep="")
	JOINED_GWAS=paste("/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_",lambda_no_dots,"/",sep="")
  QQ_PLOTS=paste("/scratch/mcostantino/descaler_gwas_plots_",lambda_no_dots,"/",sep="")
}

# path to gwas results
gwas_file <- paste0(JOINED_GWAS,pheno,"/chr1-22_",TRAIN_POP,"_",pheno,".assoc.linear")
gwas_data_load <- read.table(gwas_file, col.names=GWAS_COLS)

ci <- 0.95
nSNPs <- nrow(gwas_data_load)
print(head(gwas_data_load))

# make df w 4 cols: observed p-vals, expected p-vals, lower CI, upper CI 
# observed p-vals: sort in decreasing order, -log10 transformed - SNPs w lowest p-vals = highest value
plotdata <- data.frame(
  observed = -log10(sort(gwas_data_load$P)),
  expected = -log10(ppoints(nSNPs)), # ppoints generates sequence of probabilities given nSNPs
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))), # qbeta generates beta dist for each half of CI
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))
)

# filter SNPs w p-val > 0.01
# expected <= 2 bc -log10(0.01)=2
plotdata_sub <- plotdata %>%
  filter(expected <= 2) %>%
  sample_frac(0.01)

plotdata_sup <- plotdata %>%
  filter(expected > 2)

plotdata_small <- rbind(plotdata_sub, plotdata_sup)
print(head(plotdata_small))
rm(plotdata, plotdata_sub, plotdata_sup)

# plot the data
plot_path_vh <- paste0(QQ_PLOTS,"qqvh_",pheno,".png")
plot_path_hv <- paste0(QQ_PLOTS,"qqhv_",pheno,".png")

# qq plot
qqplot_vh <- ggplot(plotdata_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = "#00BFC4", linewidth = 1.1, direction = "vh") + 
  ggtitle(paste0(pheno, " qq plot")) +
  geom_segment(data = . %>% filter(expected == max(expected)), 
    aes(x = 0, xend = expected, y = 0, yend = expected),
    linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(
    x = expression(paste("Expected -log"[10], "(", plain(P), ")")),
    y = expression(paste("Observed -log"[10], "(", plain(P), ")"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

qqplot_hv <- ggplot(plotdata_small, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = "blue", linewidth = 1.1, direction = "hv") + 
  ggtitle(paste0(pheno," qq plot")) +
  geom_segment(data = . %>% filter(expected == max(expected)),
    aes(x = 0, xend = expected, y = 0, yend = expected),
    linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -log"[10],"(", plain(P),")")),
    y = expression(paste("Observed -log"[10],"(", plain(P),")"))) +
  theme()

ggsave(plot_path_vh, plot=qqplot_vh, width=152, height=102, units='mm')
print("plot 1 saved")

ggsave(plot_path_hv, plot=qqplot_hv, width=152, height=102, units='mm')
print("plot 2 saved")