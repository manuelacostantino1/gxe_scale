# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html

# load packages
library(dplyr)
library(ggplot2)  

#pheno=commandArgs(trailingOnly = TRUE)[1]
pheno="testosterone"
pheno_title <- gsub("_", " ", pheno)

lambda=1
print(paste("lambda=",lambda))

TRAIN_POP="whitebrit"
#JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/"
JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_male/"
GWAS_COLS <- c("CHROM","POS","ID","REF","ALT","PROVISIONAL_REF?","A1","OMITTED","A1_FREQ","TEST","OBS_CT","BETA","SE","T_STAT","P","ERRCODE")
#GWAS_PLOTS="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/gwas/"
GWAS_PLOTS="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/gwas_male/"

if(lambda==0){
	JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/"
  GWAS_PLOTS="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/plots/descaling/gwas_log/"
}

# path to gwas results
#gwas_file <- paste0(JOINED_GWAS,pheno,"/chr1-22_",TRAIN_POP,"_",pheno,".assoc.linear")
gwas_file <- paste0(JOINED_GWAS, pheno, "/chr1-22_", TRAIN_POP, "_", pheno, "_males.glm.linear")
gwas_data_load <- read.table(gwas_file, col.names=GWAS_COLS)
# ID sig & not sig snps for counts & to filter	
sig_data <- gwas_data_load %>% 
 subset(P < 0.05)
print(paste0("# sig snps (<0.05): ", dim(sig_data)[1]))
notsig_data <- gwas_data_load %>% 
  subset(P >= 0.05) %>%
  group_by(CHROM) %>% 
  sample_frac(0.1) # keep only a sample of the non-sig snps to reduce comp intensity (this filters out 10%)
gwas_data <- bind_rows(sig_data, notsig_data)

rm(gwas_data_load, sig_data, notsig_data)

# order data by creating a column with cumulative base pair position
data_cum <- gwas_data %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)
gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(bp_cum = POS + bp_add)

# get centre position of each chromosome
axis_set <- gwas_data %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))
# set the limit of the y-axis - avoid cutting off highly significant SNPs	
ylim <- gwas_data %>% 
  filter(P == min(P)) %>% 
  mutate(ylim = abs(floor(log10(P))) + 2) %>% 
  pull(ylim)
test1 <- gwas_data %>%
  filter(P == min(P))
test2 <- gwas_data %>%
  filter(P == min(P)) %>%
  mutate(ylim = log10(P))

# set Bonferroni threshold
sig <- 5e-8

# plot the data
plot_path <- paste0(GWAS_PLOTS,"manhattan_",pheno,".png")

cmplot_colors <- c("deepskyblue", "orange", "mediumseagreen", "orchid", "gold",
                   "blue", "green", "red", "purple", "brown")

# Extend the palette to match number of chromosomes (recycle if needed)
num_chr <- length(unique(gwas_data$CHROM))
chr_palette <- rep(cmplot_colors, length.out = num_chr)

# Plot
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P), color = as.factor(CHROM))) +
  geom_hline(yintercept = -log10(5e-8), color = "grey40", linetype = "dashed") +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = chr_palette) +
  scale_x_continuous(labels = axis_set$CHROM, breaks = axis_set$center) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = paste0(pheno, " GWAS")
  ) +
  theme_minimal(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
print("done")
ggsave(plot_path, plot=manhplot, width=152, height=102, units='mm')

print("plot saved")	
rm(ylim, axis_set, data_cum, gwas_data)