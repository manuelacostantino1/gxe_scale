rm(list = ls())
set.seed(34)

library(ggplot2)
library(dplyr)
library(tidyr)

# Quantile normalization 
quantnorm <- function(Y) {
  obs <- which(!is.na(Y))
  ranks <- rank(Y[obs])
  quantiles <- qnorm((ranks + 1) / (length(ranks) + 2))
  if (!all.equal(ranks, rank(quantiles))) stop("Bad rank problem")
  Y[obs] <- quantiles
  scale(Y)
}

# Parameters
niter <- 500
N <- 1e5
P <- 25       
lamseq <- seq(-1, 2, length.out = 25)

# Data frame to store data
scenarios <- data.frame(
  scenario = c(
    "Linear scale - GxE",
    "Exponential scale - GxE",
    "Inverse-logit scale - GxE",
    "Linear scale - no GxE",
    "Exponential scale - no GxE",
    "Inverse-logit scale - no GxE"
  ),
  truelam = c(1, 0, 100, 1, 0, 100),
  gxe_var = c(0.1, 0.1, 0.1, 0, 0, 0),
  stringsAsFactors = FALSE
)

# Simulation function
run_simulation <- function(niter, N, P, lamseq, scenarios, savefile) {
  if (file.exists(savefile)) {
    message("Loading data from: ", savefile)
    load(savefile)
  } else {
    message("Running simulation because no file found...")
    
    # Storage
    pvs <- array(dim = c(nrow(scenarios), length(lamseq), P, niter))
    pvs_qn <- array(dim = c(nrow(scenarios), P, niter))
    pvs_untrans <- array(dim = c(nrow(scenarios), P, niter))
    
    for (it in 1:niter) {
      cat("Iteration", it, "\n")
      for (i in seq_len(nrow(scenarios))) {
        truelam <- scenarios$truelam[i]
        gxe_var <- scenarios$gxe_var[i]
        
        # Simulate G and E
        G <- matrix(rnorm(N * P), N, P)
        e <- scale(rbinom(N, 1, 0.3))
        
        # Variance components
        g_var <- 0.3
        e_var <- 0.2
        resid_var <- 1 - (g_var + gxe_var + e_var)
        
        # Phenotype
        y <- scale(rowSums(G)) * sqrt(g_var) +
          scale((G * (as.numeric(e) %o% rep(1, P))) %*% rnorm(P)) * sqrt(gxe_var) +
          e * sqrt(e_var) +
          rnorm(N) * sqrt(resid_var)
        
        y <- y + (0.1 * sd(y) - min(y))  # shift positive
        
        # Apply true transformation
        if (truelam == 0) {
          y <- exp(y)
        } else if (truelam == 100) {
          y <- exp(y) / (1 + exp(y))
        } else {
          y <- y^(1 / truelam)
        }
        
        # Candidate lambdas
        for (l in seq_along(lamseq)) {
          lam <- lamseq[l]
          yl <- if (lam == 0) scale(log(y)) else scale((y^lam - 1) / lam)
          pvs[i, l, , it] <- apply(G, 2, function(g) summary(lm(yl ~ g * e))$coef["g:e", 4])
        }
        
        # Quantile-normalized p-values
        yqn <- quantnorm(y)
        pvs_qn[i, , it] <- apply(G, 2, function(g) summary(lm(yqn ~ g * e))$coef["g:e", 4])
        
        # Untransformed lambda=1 p-values
        lambda_1_index <- which.min(abs(lamseq - 1))
        pvs_untrans[i, , it] <- pvs[i, lambda_1_index, , it]
      }
    }
    
    # Dimnames
    dimnames(pvs) <- list(
      scenario = scenarios$scenario,
      lambda = as.character(lamseq),
      snp = paste0("SNP", 1:P),
      iter = paste0("iter", 1:niter)
    )
    dimnames(pvs_qn) <- list(
      scenario = scenarios$scenario,
      snp = paste0("SNP", 1:P),
      iter = paste0("iter", 1:niter)
    )
    dimnames(pvs_untrans) <- list(
      scenario = scenarios$scenario,
      snp = paste0("SNP", 1:P),
      iter = paste0("iter", 1:niter)
    )
    
    save(pvs, pvs_qn, pvs_untrans, scenarios, lamseq, P, niter, file = savefile)
    message("Saved simulation results to: ", savefile)
  }
  
  list(pvs = pvs, pvs_qn = pvs_qn, pvs_untrans = pvs_untrans)
}

# Run
savefile <- "/gpfs/data/ukb-share/dahl/manuela/descaler_sims/sims_pvs_results.RData"
res <- run_simulation(niter, N, P, lamseq, scenarios, savefile)

# Convert to long format for plotting
pvs_long <- as.data.frame.table(res$pvs, responseName = "p", stringsAsFactors = FALSE) %>%
  mutate(
    lambda = as.numeric(lambda),
    iter = as.integer(gsub("iter", "", iter)),
    log10p = -log10(p)
  ) %>%
  mutate(log10p = ifelse(is.infinite(log10p), max(log10p[is.finite(log10p)], na.rm = TRUE), log10p))

# Summarize mean across SNPs
pvs_snp_summary <- pvs_long %>%
  group_by(scenario, lambda, iter) %>%
  summarize(mean_log10p = mean(log10p), .groups = "drop")

# Reorder panels for plot
scenario_order <- c(
  "Linear scale - no GxE",
  "Exponential scale - no GxE",
  "Inverse-logit scale - no GxE",
  "Linear scale - GxE",
  "Exponential scale - GxE",
  "Inverse-logit scale - GxE"
)

pvs_snp_summary$scenario <- factor(pvs_snp_summary$scenario, levels = scenario_order)

# Plot
#png("~/Desktop/simulation_line.png", width = 8, height = 5, units = "in", res = 400)
#ggplot(pvs_snp_summary, aes(x = lambda, y = mean_log10p, group = iter)) +
#  geom_line(aes(color = "Each simulation"), alpha = 0.4) +
#  geom_hline(aes(yintercept = -log10(0.05 / P), color = "Significance threshold"),
#             linetype = "dashed", size = 0.8) +
#  scale_color_manual(name = "Legend",
#                     values = c("Each simulation" = "#00CCFF",
#                                "Significance threshold" = "black")) +
#  facet_wrap(~scenario, nrow = 2, scales = "free_y") +
#  labs(x = expression(paste("Box-Cox ", lambda)),
#       y = expression(-log[10](p))) +
#  theme_minimal(base_size = 10) +
#  theme(
#    axis.text = element_text(size = 8),
#    axis.title = element_text(size = 10),
#    legend.text = element_text(size = 8),
#    legend.title = element_text(size = 9),
#    strip.text = element_text(size = 9, face = "bold"),
#    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
#  )+
#  guides(color = guide_legend(override.aes = list(linetype = 1)))
#dev.off()




# Calculate FPR in QN and untransformed datasets
threshold <- 0.05
dim(res$pvs_untrans)
fpr_untrans <- sapply(seq_len(nrow(scenarios)), function(i) {
  pvals <- res$pvs_untrans[i, , ]   
  mean(pvals < threshold)
})

fpr_qn <- sapply(seq_len(nrow(scenarios)), function(i) {
  pvals <- res$pvs_qn[i, , ]   
  mean(pvals < threshold)
})

fpr_df <- data.frame(
  scenario = scenarios$scenario,
  FPR_untrans = fpr_untrans,
  FPR_qn = fpr_qn
)

fpr_df

write.table(fpr_df, file="/gpfs/data/ukb-share/dahl/manuela/descaler_sims/rint_fpr.txt", row.names=FALSE, quote=FALSE, sep="\t")

