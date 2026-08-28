# vQTL sensitivity simulation
# Uses helper functions in sims/sim_helpers.R

source("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/descaler/sims/sim_helpers.R")

library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
threshold <- 0.05

# Allow quick pilot runs by setting pilot=TRUE
pilot <- FALSE
if (pilot) {
  niter <- 5
  N <- 1000
  P <- 3
} else {
  niter <- 500
  N <- 1e5
  P <- 25
}

# Run one vqtl_var value across four methods: untransformed, RINT, log, and sqrt (lambda=0.5)
# observe_exp: if TRUE, exponentiate y before testing (panel 2)
run_vqtl <- function(vqtl_var, observe_exp = FALSE) {
  total_tests <- niter * P
  count_untrans <- 0
  count_qn <- 0
  count_log <- 0
  count_sqrt <- 0

  for (it in 1:niter) {
    G <- matrix(rnorm(N * P), N, P)
    e <- scale(rbinom(N, 1, 0.3))

    y <- generate_phenotype(
      G = G,
      e = e,
      g_var = 0.3,
      e_var = 0.2,
      gxe_var = 0,
      vqtl_var = vqtl_var
    )

    y_obs <- if (observe_exp) exp(y) else y

    # Untransformed
    p_untrans <- apply(G, 2, function(g) summary(lm(y_obs ~ g * e))$coef["g:e", 4])
    count_untrans <- count_untrans + sum(p_untrans < threshold, na.rm = TRUE)

    # RINT (quantile normalisation)
    y_qn <- transform_qn(y_obs)
    p_qn <- apply(G, 2, function(g) summary(lm(y_qn ~ g * e))$coef["g:e", 4])
    count_qn <- count_qn + sum(p_qn < threshold, na.rm = TRUE)

    # Log (lambda = 0)
    y_log <- transform_boxcox(y_obs, lambda = 0)
    p_log <- apply(G, 2, function(g) summary(lm(y_log ~ g * e))$coef["g:e", 4])
    count_log <- count_log + sum(p_log < threshold, na.rm = TRUE)

    # Sqrt (lambda = 0.5)
    y_sqrt <- transform_boxcox(y_obs, lambda = 0.5)
    p_sqrt <- apply(G, 2, function(g) summary(lm(y_sqrt ~ g * e))$coef["g:e", 4])
    count_sqrt <- count_sqrt + sum(p_sqrt < threshold, na.rm = TRUE)
  }

  data.frame(
    vqtl_var = vqtl_var,
    FPR_untrans = count_untrans / total_tests,
    FPR_qn = count_qn / total_tests,
    FPR_log = count_log / total_tests,
    FPR_sqrt = count_sqrt / total_tests
  )
}

vqtl_var_seq <- c(0, 0.05, 0.1, 0.15, 0.2)

outdir <- "/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/intermediate_files"
outdir <- path.expand(outdir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

rdata_path <- file.path(outdir, "vqtl_sweep_results.RData")

if (file.exists(rdata_path)) {
  cat("Loading existing results from", rdata_path, "\n")
  load(rdata_path)
} else {
  results_same <- bind_rows(lapply(vqtl_var_seq, run_vqtl, observe_exp = FALSE))
  results_exp  <- bind_rows(lapply(vqtl_var_seq, run_vqtl, observe_exp = TRUE))
  results_same$panel <- "Observe on generated scale"
  results_exp$panel  <- "Observe on exp(generated scale)"
  results <- bind_rows(results_same, results_exp)

  save(results, results_same, results_exp, vqtl_var_seq, threshold,
       file = rdata_path)
  write.table(results, file = file.path(outdir, "vqtl_sweep_summary.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Plot: FPR vs vQTL strength for three methods, two panels
plot_long <- results %>%
  pivot_longer(
    cols = starts_with("FPR_"),
    names_to = "method",
    values_to = "FPR"
  ) %>%
  mutate(
    method = case_when(
      method == "FPR_untrans" ~ "Untransformed",
      method == "FPR_qn" ~ "RINT",
      method == "FPR_log" ~ "Log",
      method == "FPR_sqrt" ~ "Sqrt~(lambda==0.5)",
      TRUE ~ method
    ),
    panel = factor(panel, levels = c("Observe on generated scale", "Observe on exp(generated scale)"))
  )

p <- ggplot(plot_long, aes(x = vqtl_var, y = FPR, color = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~panel) +
  scale_color_manual(
    values = c(
      "Untransformed"      = "#0072B2",
      "RINT"               = "#009E73",
      "Log"                = "#D55E00",
      "Sqrt~(lambda==0.5)" = "#CC79A7"
    ),
    labels = function(x) parse(text = x)
  ) +
  labs(
    x = "Variance explained by vQTL",
    y = "False Positive Rate (GxE detection)",
    color = "Method"
  ) +
  theme_minimal(base_size = 11)

png(file.path(outdir, "vqtl_sweep_fpr.png"), width = 12, height = 4.5, units = "in", res = 200)
print(p)
dev.off()

cat("Saved vQTL sweep results to:", outdir, "\n")
