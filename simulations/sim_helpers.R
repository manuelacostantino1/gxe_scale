# Helper functions for sim+qn.R moved here so they can be sourced independently

# Quantile normalization
quantnorm <- function(Y) {
  obs <- which(!is.na(Y))
  ranks <- rank(Y[obs])
  quantiles <- qnorm((ranks + 1) / (length(ranks) + 2))
  if (!all.equal(ranks, rank(quantiles))) stop("Bad rank problem")
  Y[obs] <- quantiles
  scale(Y)
}

# Transformation helpers
transform_identity <- identity

transform_boxcox <- function(y, lambda) {
  if (lambda == 0) {
    scale(log(y))
  } else {
    scale((y^lambda - 1) / lambda)
  }
}

transform_log_offset <- function(y, c_val) {
  scale(log(y + c_val))
}

transform_qn <- quantnorm

# Phenotype generator (matches original construction; supports optional vQTL residual heteroskedasticity)
# vQTL variance modifier is rowMeans(G)
generate_phenotype <- function(G, e, g_var = 0.3, e_var = 0.2, gxe_var = 0, vqtl_var = 0, truelam = 1) {
  N <- nrow(G)
  P <- ncol(G)
  resid_var <- 1 - (g_var + gxe_var + e_var + vqtl_var)

  e_latent <- rnorm(N)

  y <- scale(rowSums(G)) * sqrt(g_var) +
    scale((G * (as.numeric(e) %o% rep(1, P))) %*% rnorm(P)) * sqrt(gxe_var) +
    scale(rowMeans(G) * e_latent) * sqrt(vqtl_var) +
    e * sqrt(e_var) +
    rnorm(N) * sqrt(resid_var)

  y <- y + (0.1 * sd(y) - min(y))  # shift positive

  if (truelam == 0) {
    y <- exp(y)
  } else if (truelam == 100) {
    y <- exp(y) / (1 + exp(y))
  } else {
    y <- y^(1 / truelam)
  }

  y
}

# Extract p-values for GxE across Box-Cox lambdas, QN and untransformed
extract_pvalues <- function(G, e, y, lamseq) {
  P <- ncol(G)
  # matrix: lambda x SNP
  pvs <- matrix(NA, nrow = length(lamseq), ncol = P)
  for (l in seq_along(lamseq)) {
    lam <- lamseq[l]
    yl <- if (lam == 0) scale(log(y)) else scale((y^lam - 1) / lam)
    pvs[l, ] <- apply(G, 2, function(g) summary(lm(yl ~ g * e))$coef["g:e", 4])
  }

  # quantile-normalized p-values
  yqn <- transform_qn(y)
  pvs_qn <- apply(G, 2, function(g) summary(lm(yqn ~ g * e))$coef["g:e", 4])

  # untransformed (lambda ~= 1)
  lambda_1_index <- which.min(abs(lamseq - 1))
  pvs_untrans <- pvs[lambda_1_index, ]

  list(pvs = pvs, pvs_qn = pvs_qn, pvs_untrans = pvs_untrans)
}


