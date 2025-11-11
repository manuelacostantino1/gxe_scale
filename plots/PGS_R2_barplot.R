library(ggplot2)
library(tidyr)
library(dplyr)
pheno <- "height"

# THE SEXES ARE ALSO FLIPPED BC THERE IS AN ERRROR UPSTREAM THAT NEEDS TO BE FIXED

# Open files 
data_female <- read.table(paste0("/Users/manuelacostantino/Desktop/descaler_analysis/log_vs_lin/log_vs_lin_",pheno,"_white_euro_male.txt"))
data_male <- read.table(paste0("/Users/manuelacostantino/Desktop/descaler_analysis/log_vs_lin/log_vs_lin_",pheno,"_white_euro_female.txt"))
data_all <- read.table(paste0("/Users/manuelacostantino/Desktop/descaler_analysis/log_vs_lin/log_vs_lin_",pheno,"_white_euro_all.txt"))
head(data_male)
data_all$group <- "All"
data_male$group <- "Male"
data_female$group <- "Female"

# Use the first row from each group (assumed summary row)
df <- rbind(data_all[1,], data_male[1,], data_female[1,])

# Reshape and relabel for plotting
df_long <- df %>%
  select(group, r20, r21) %>%
  rename(`Log Scale` = r20, `Default Scale` = r21) %>%
  pivot_longer(cols = c(`Log Scale`, `Default Scale`), names_to = "Scale", values_to = "Value")

# Fix legend label order
df_long$Scale <- factor(df_long$Scale, levels = c("Log Scale", "Default Scale"))

# Plot
png(paste0("~/Desktop/",pheno,"_r2.png"), width = 3, height = 3,units="in",res=300)
ggplot(df_long, aes(x = group, y = Value, fill = Scale)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6,
    color = "black",
    size = 0.3   
  ) +
  labs(
    x = NULL,
    y = expression(R^2),
    fill = "Scale"
  ) +
  scale_fill_manual(
    values = c("Log Scale" = "mediumblue", "Default Scale" = "lightcyan2"),
    labels = c("Log", "Default")
  ) +
  theme_minimal(base_size = 8) +  
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),  
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9),       
    axis.text.x = element_text(size = 7),     
    axis.text.y = element_text(size = 7),    
    legend.text = element_text(size = 7),      
    legend.title = element_text(size = 8),     
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 9) 
  )
dev.off()
