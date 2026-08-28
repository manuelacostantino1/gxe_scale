# Phase 7: vQTL detection power/FPR across Box-Cox transformations
# Brown-Forsythe test on discrete genotypes G ~ Binomial(2, maf)
# Mirrors sim+qn.R but tests variance heterogeneity (vQTL), not GxE

source("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/descaler/sims/sim_helpers.R")

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

pilot <- TRUE
if (pilot) {
  niter <- 5
  N     <- 1000
  P     <- 3
} else {
  niter <- 500
  N     <- 1e5
  P     <- 25
}

lamseq     <- seq(-1, 2, length.out = 25)
maf        <- 0.3
vqtl_alpha <- 1.0  # heteroskedasticity strength for vQTL scenarios

scenarios <- data.frame(
  scenario = c(
    "Linear scale - vQTL",
    "Exponential scale - vQTL",
    "Inverse-logit scale - vQTL",
    "Linear scale - no vQTL",
    "Exponential scale - no vQTL",
    "Inverse-logit scale - no vQTL"
  ),
  truelam  = c(1, 0, 100, 1, 0, 100),
  has_vqtl = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Efficient Brown-Forsythe test (Levene's test with group medians)
# Returns p-value for variance heterogeneity across genotype groups {0, 1, 2}
bf_pvalue <- function(y, g) {
  gvals <- sort(unique(g))
  if (length(gvals) < 2) return(NA_real_)
  n    <- length(y)
  k    <- length(gvals)
  nj   <- as.integer(table(factor(g, levels = gvals)))
  meds <- tapply(y, g, median)
  z    <- abs(y - meds[as.character(g)])
  zbar_j <- tapply(z, g, mean)
  zbar   <- mean(z)
  SS_b <- sum(nj * (zbar_j - zbar)^2)
  SS_w <- sum(tapply(z, g, function(zi) sum((zi - mean(zi))^2)))
  if (SS_w == 0) return(NA_real_)
  F_stat <- (SS_b / (k - 1)) / (SS_w / (n - k))
  pf(F_stat, df1 = k - 1, df2 = n - k, lower.tail = FALSE)
}

outdir   <- "/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/intermediate_files"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
savefile <- file.path(outdir, "phase7_vqtl_detection_results.RData")

if (file.exists(savefile)) {
  message("Loading saved results from ", savefile)
  load(savefile)
} else {
  message("Running Phase 7 simulation...")

  # [scenario x lambda x SNP x iter]
  pvs <- array(dim = c(nrow(scenarios), length(lamseq), P, niter))

  for (it in seq_len(niter)) {
    cat("Iteration", it, "\n")
    for (i in seq_len(nrow(scenarios))) {
      truelam  <- scenarios$truelam[i]
      has_vqtl <- scenarios$has_vqtl[i]

      # Discrete genotypes: mimic SNP dosage
      G <- matrix(rbinom(N * P, 2, maf), N, P)

      g_var     <- 0.3
      resid_var <- 1 - g_var

      y_mean <- as.numeric(scale(rowSums(G))) * sqrt(g_var)

      if (has_vqtl) {
        # Polygenic vQTL: residual variance scales with rowMeans(G)
        vqtl_score <- abs(as.numeric(scale(rowMeans(G))))  # U-shaped: extremes (G=0,G=2) get more variance than middle (G=1)
        sd_factor  <- exp(vqtl_alpha * vqtl_score / 2)
        eps        <- rnorm(N) * sd_factor
        eps        <- eps / sd(eps) * sqrt(resid_var)
      } else {
        eps <- rnorm(N) * sqrt(resid_var)
      }

      y <- y_mean + eps

      if (truelam == 100) {
        y <- exp(y) / (1 + exp(y))   # plogis on latent y ~ N(0,1); no shift needed
      } else {
        y <- y + (0.1 * sd(y) - min(y))  # shift positive
        if (truelam == 0) {
          y <- exp(y)
        } else {
          y <- y^(1 / truelam)
        }
      }

      for (l in seq_along(lamseq)) {
        lam <- lamseq[l]
        yl  <- if (lam == 0) as.numeric(scale(log(y))) else as.numeric(scale((y^lam - 1) / lam))
        pvs[i, l, , it] <- vapply(seq_len(P), function(j) bf_pvalue(yl, G[, j]), numeric(1))
      }
    }
  }

  dimnames(pvs) <- list(
    scenario = scenarios$scenario,
    lambda   = as.character(lamseq),
    snp      = paste0("SNP", seq_len(P)),
    iter     = paste0("iter", seq_len(niter))
  )

  save(pvs, scenarios, lamseq, P, niter, file = savefile)
  message("Saved simulation results to: ", savefile)
}

# Long format
pvs_long <- as.data.frame.table(pvs, responseName = "p", stringsAsFactors = FALSE) %>%
  mutate(
    lambda = as.numeric(lambda),
    iter   = as.integer(gsub("iter", "", iter)),
    log10p = -log10(p)
  ) %>%
  mutate(log10p = ifelse(is.infinite(log10p), max(log10p[is.finite(log10p)], na.rm = TRUE), log10p))

# Mean across SNPs per scenario/lambda/iter
pvs_snp_summary <- pvs_long %>%
  group_by(scenario, lambda, iter) %>%
  summarize(mean_log10p = mean(log10p, na.rm = TRUE), .groups = "drop")

scenario_order <- c(
  "Linear scale - no vQTL",
  "Exponential scale - no vQTL",
  "Inverse-logit scale - no vQTL",
  "Linear scale - vQTL",
  "Exponential scale - vQTL",
  "Inverse-logit scale - vQTL"
)
pvs_snp_summary$scenario <- factor(pvs_snp_summary$scenario, levels = scenario_order)

png(file.path(outdir, "phase7_vqtl_detection.png"), width = 8, height = 5, units = "in", res = 400)
ggplot(pvs_snp_summary, aes(x = lambda, y = mean_log10p, group = iter)) +
  geom_line(aes(color = "Each simulation"), alpha = 0.4) +
  geom_hline(aes(yintercept = -log10(0.05 / P), color = "Significance threshold"),
             linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(name = "Legend",
                     values = c("Each simulation" = "#00CCFF",
                                "Significance threshold" = "black")) +
  facet_wrap(~scenario, nrow = 2, scales = "free_y") +
  labs(
    x = expression(paste("Box-Cox ", lambda)),
    y = expression(-log[10](p))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 10),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    strip.text   = element_text(size = 9, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1)))
dev.off()

cat("Phase 7 complete. Outputs saved to:", outdir, "\n")
