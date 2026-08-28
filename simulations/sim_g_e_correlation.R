# Correlated G-E simulation
source("/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/descaler/sims/sim_helpers.R")

library(ggplot2)
library(dplyr)
library(tidyr)

threshold <- 0.05
rho_seq   <- c(0, 0.1, 0.3, 0.5)

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

# Inverse of the three observation-scale transforms:
#   NA  -> identity (linear scale)
#   0   -> exp()    (log scale)
#   0.5 -> Box-Cox inverse with lambda=0.5 (sqrt scale)
bc_inv <- function(y, lambda) {
  if (is.na(lambda)) y
  else if (lambda == 0) exp(y)
  else (lambda * y + 1)^(1 / lambda)
}

obs_scales <- list(
  list(label = "Linear scale", lambda = NA),
  list(label = "Log scale",    lambda = 0),
  list(label = "Sqrt scale",   lambda = 0.5)
)

make_correlated_environment <- function(G, rho) {
  g_signal <- as.numeric(scale(rowMeans(G)))
  e_noise  <- as.numeric(scale(rbinom(nrow(G), 1, 0.3)))
  as.numeric(scale(rho * g_signal + sqrt(1 - rho^2) * e_noise))
}

fit_gxe_stats <- function(response, G, e) {
  P <- ncol(G)
  p_values <- numeric(P)
  betas    <- numeric(P)
  for (snp in seq_len(P)) {
    df  <- data.frame(y = response, g = G[, snp], e = e)
    fit <- lm(y ~ g * e, data = df)
    ct  <- summary(fit)$coef
    betas[snp]    <- ct["g:e", "Estimate"]
    p_values[snp] <- ct["g:e", "Pr(>|t|)"]
  }
  list(p_values = p_values, betas = betas)
}

run_scenario <- function(gxe_var, scenario_name, obs_label, obs_lambda) {
  method_names <- c("untrans", "rint", "log", "sqrt")
  results <- vector("list", length(rho_seq))

  for (rho_idx in seq_along(rho_seq)) {
    rho <- rho_seq[rho_idx]
    cat(sprintf("  %s | %s | rho=%.1f\n", scenario_name, obs_label, rho))

    method_counts <- setNames(rep(0L, 4L), method_names)
    method_betas  <- setNames(replicate(4, numeric(0), simplify = FALSE), method_names)

    for (it in seq_len(niter)) {
      G      <- matrix(rnorm(N * P), N, P)
      e      <- make_correlated_environment(G, rho)
      y_base <- as.numeric(generate_phenotype(G, e, g_var = 0.3, e_var = 0.2,
                                              gxe_var = gxe_var, truelam = 1))
      # exp(y_base) would be exp(N(4.5,1)) ~ lognormal with mean≈90, breaking OLS;
      # center first so exp(N(0,1)) gives a manageable lognormal (mean≈1.65)
      y_in  <- if (!is.na(obs_lambda) && obs_lambda == 0) as.numeric(scale(y_base)) else y_base
      y_obs <- bc_inv(y_in, obs_lambda)

      ms <- list(
        untrans = fit_gxe_stats(y_obs, G, e),
        rint    = fit_gxe_stats(as.numeric(transform_qn(y_obs)), G, e),
        log     = fit_gxe_stats(as.numeric(transform_boxcox(y_obs, 0)), G, e),
        sqrt    = fit_gxe_stats(as.numeric(transform_boxcox(y_obs, 0.5)), G, e)
      )

      for (mn in method_names) {
        method_counts[[mn]] <- method_counts[[mn]] + sum(ms[[mn]]$p_values < threshold, na.rm = TRUE)
        method_betas[[mn]]  <- c(method_betas[[mn]], ms[[mn]]$betas)
      }
    }

    total_tests <- niter * P
    results[[rho_idx]] <- data.frame(
      scenario  = scenario_name,
      obs_scale = obs_label,
      rho       = rho,
      method    = c("Untransformed", "RINT", "Log", "Sqrt~(lambda==0.5)"),
      FPR       = c(method_counts[["untrans"]],
                    method_counts[["rint"]],
                    method_counts[["log"]],
                    method_counts[["sqrt"]]) / total_tests,
      beta_mean = c(mean(method_betas[["untrans"]]),
                    mean(method_betas[["rint"]]),
                    mean(method_betas[["log"]]),
                    mean(method_betas[["sqrt"]])),
      stringsAsFactors = FALSE
    )
  }

  bind_rows(results)
}

outdir <- "/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/intermediate_files"
outdir <- path.expand(outdir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

rdata_path <- file.path(outdir, "ge_corr_sweep_results.RData")

if (file.exists(rdata_path)) {
  cat("Loading existing results from", rdata_path, "\n")
  load(rdata_path)
} else {
  # Run all 6 combinations: 3 obs_scales x 2 scenarios
  all_results <- lapply(obs_scales, function(os) {
    bind_rows(
      run_scenario(0,   "null",   os$label, os$lambda),
      run_scenario(0.1, "signal", os$label, os$lambda)
    )
  })
  results <- bind_rows(all_results)

  # beta_bias: deviation from rho=0 estimate within each scenario/scale/method
  results <- results %>%
    group_by(scenario, obs_scale, method) %>%
    mutate(beta_bias = beta_mean - beta_mean[rho == 0]) %>%
    ungroup()

  save(results, rho_seq, threshold, file = rdata_path)
  write.table(results, file = file.path(outdir, "ge_corr_sweep_summary.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Build 3-column plot data
null_res   <- filter(results, scenario == "null")
signal_res <- filter(results, scenario == "signal")

plot_data <- bind_rows(
  null_res   %>% transmute(obs_scale, rho, method,
                            metric = "False Positive Rate",
                            value  = FPR),
  signal_res %>% transmute(obs_scale, rho, method,
                            metric = "Coefficient deviation\n(GxE present)",
                            value  = beta_bias),
  null_res   %>% transmute(obs_scale, rho, method,
                            metric = "Coefficient deviation\n(GxE absent)",
                            value  = beta_mean)
) %>%
  mutate(
    obs_scale = factor(obs_scale, levels = sapply(obs_scales, `[[`, "label")),
    metric    = factor(metric, levels = c("False Positive Rate",
                                          "Coefficient deviation\n(GxE present)",
                                          "Coefficient deviation\n(GxE absent)"))
  )

library(grid)
library(gtable)

color_vals <- c(
  "Untransformed"      = "#0072B2",
  "RINT"               = "#009E73",
  "Log"                = "#D55E00",
  "Sqrt~(lambda==0.5)" = "#CC79A7"
)

# Build one ggplot per obs_scale row; facet_wrap over metric gives free y per cell
make_row <- function(obs_lbl, show_xlabel, show_legend) {
  ggplot(filter(plot_data, obs_scale == obs_lbl),
         aes(x = rho, y = value, color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = color_vals, labels = function(x) parse(text = x)) +
    labs(x = if (show_xlabel) expression(rho ~ "(G-E correlation)") else NULL,
         y = NULL, color = "Method") +
    theme_minimal(base_size = 11) +
    theme(
      strip.text      = element_text(size = 9, face = "bold"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position  = if (show_legend) "top" else "none",
      legend.direction = "horizontal"
    )
}

# Append a right-side strip label (obs_scale row name) to a gtable — no background
add_right_strip <- function(g, label) {
  pr <- g$layout[grepl("^panel", g$layout$name), ]
  g  <- gtable_add_cols(g, unit(0.6, "cm"), pos = ncol(g))
  gtable_add_grob(g,
    textGrob(label, rot = -90, gp = gpar(fontsize = 9, fontface = "bold")),
    t = min(pr$t), b = max(pr$b),
    l = ncol(g),   r = ncol(g),
    name = paste0("right-strip-", label))
}

# First row carries the legend at top; others suppress it
g1 <- add_right_strip(ggplotGrob(make_row("Linear scale", FALSE, TRUE)),  "Linear scale")
g2 <- add_right_strip(ggplotGrob(make_row("Log scale",    FALSE, FALSE)), "Log scale")
g3 <- add_right_strip(ggplotGrob(make_row("Sqrt scale",   TRUE,  FALSE)), "Sqrt scale")

# Align column widths so panels line up vertically
maxw <- unit.pmax(g1$widths, g2$widths, g3$widths)
g1$widths <- g2$widths <- g3$widths <- maxw

png(file.path(outdir, "ge_corr_sweep.png"), width = 10, height = 8, units = "in", res = 200)
grid.newpage()
grid.draw(rbind(g1, g2, g3, size = "first"))
dev.off()

cat("Saved G-E correlation sweep results to:", outdir, "\n")
