# Phase 6: Generate from inverse Box-Cox and test with log(y + c)

source("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/descaler/sims/sim_helpers.R")

library(ggplot2)
library(dplyr)
library(tidyr)

# Parameters
c_values   <- c(0, 1, 10, 100)
gen_lambdas <- c(0, 0.25, 0.5, 0.75, 1)

# Pilot vs full
pilot <- FALSE
if (pilot) {
  niter <- 10
  N <- 1000
  P <- 3
} else {
  niter <- 500
  N <- 1e5
  P <- 25
}

bc_inv <- function(y, lambda) {
  if (lambda == 0) exp(y) else (lambda * y + 1)^(1 / lambda)
}

run_phase6 <- function(gxe_var, scenario_name) {
  res_all <- list()

  make_positive <- function(x) x - min(x, na.rm = TRUE) + 1e-6

  for (gen_l in gen_lambdas) {
    agg_array <- array(0, dim = c(niter, length(c_values)))

    for (it in seq_len(niter)) {
      G <- matrix(rnorm(N * P), N, P)
      e <- scale(rbinom(N, 1, 0.3))

      y_base <- generate_phenotype(
        G = G, e = e,
        g_var = 0.3, e_var = 0.2,
        gxe_var = gxe_var,
        truelam = 1
      )

      # Generate observed y as inverse Box-Cox of y_base
      # bc_inv of positive y_base is always positive — no make_positive needed
      y_obs <- bc_inv(as.numeric(y_base), gen_l)

      for (ci in seq_along(c_values)) {
        c_val <- c_values[ci]
        y_test_raw <- log(y_obs + c_val)
        y_test <- as.numeric(scale(y_test_raw))
        pvs <- apply(G, 2, function(g) summary(lm(y_test ~ g * e))$coef["g:e", 4])
        mlogp <- median(-log10(pvs), na.rm = TRUE)
        agg_array[it, ci] <- mlogp
      }
    }

    mean_mlogp <- apply(agg_array, 2, mean, na.rm = TRUE)
    df <- data.frame(
      c        = c_values,
      mean_mlogp = mean_mlogp,
      gen_lambda = gen_l,
      scenario = scenario_name
    )
    res_all[[as.character(gen_l)]] <- df
  }

  bind_rows(res_all)
}

results_null <- run_phase6(gxe_var = 0, scenario_name = "No GxE")

outdir <- path.expand("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/intermediate_files")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

save(results_null, file = file.path(outdir, "phase6_inverse_log_boxcox_results.RData"))
write.table(results_null, file.path(outdir, "phase6_inverse_null_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_df <- results_null %>%
  mutate(gen_lambda = factor(gen_lambda, levels = as.character(gen_lambdas)))

p <- ggplot(plot_df, aes(x = c, y = mean_mlogp, color = gen_lambda, group = gen_lambda)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_color_viridis_d(option = "plasma", begin = 0.15, end = 0.9) +
  labs(
    x     = "Test log offset (c)",
    y     = "median -log10(p) (g:e)",
    color = "Generation Box-Cox lambda"
  ) +
  theme_minimal(base_size = 11)

png(file.path(outdir, "phase6_log_inverse_lambdas_by_c.png"),
    width = 6, height = 4.5, units = "in", res = 200)
print(p)
dev.off()

cat("Saved Phase 6 results to:", outdir, "\n")
