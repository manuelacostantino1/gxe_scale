# calculate R2
r2fxn <- function(x,y,X){
 X <- as.matrix(X)
 y <- resid( lm( y ~ X ,na.action = na.exclude) )
 x <- resid( lm( x ~ X ,na.action = na.exclude) )
 summary( lm( y ~ x ) )
}

r2fxn_spearman <- function(x,y,X){
 X <- as.matrix(X)
 y <- resid( lm( y ~ X ,na.action = na.exclude) )
 x <- resid( lm( x ~ X ,na.action = na.exclude) )
 (cor(x, y, method = "spearman", use="complete.obs"))^2
}

r2_diff_spearman <- function(x0,x1,y,B=1e4,X){
 r21 <- r2fxn_spearman(x1,y,X)
 r20 <- r2fxn_spearman(x0,y,X)
 delta <- r21-r20
 delta_boots <- numeric()
 for( b in 1:B ){
  is <- sample( 1:length(x0), replace=T )
  x0b <- x0[is]
  x1b <- x1[is]
  yb <- y[is]
  Xb <- X[is,]
  delta_boots[b] <- (r2fxn_spearman(x0b,yb,Xb)-r2fxn_spearman(x1b,yb,Xb)) / r2fxn_spearman(x1b,yb,Xb)
 }
 pv_oneside <- (1 + sum(delta_boots <= 0))/ (1+B)
 pval    <- 2*min(pv_oneside,1 - pv_oneside)
 list( r20=r20, r21=r21, ci_95=quantile(delta_boots,c(0.025,.975)), pval=pval)
}

# bootstrap R2
r2_boot <- function(x,y,B=1e4,X){
 ## x=pgs, y=pheno, X=covars
 lm_out <- r2fxn(x,y,X)
 coef <- ifelse("x" %in% rownames(lm_out$coefficients), lm_out$coefficients["x", "Estimate"], NA)
 if (is.na(coef)) {cat("Error: Coefficient 'x' not found in lm_out, pgs might be 0s\n")}
 r2 <- lm_out$r.squared
 r2_boots <- NULL
 for( b in 1:B ){
  is <- sample( 1:length(x), replace=T )
  x1 <- x[is]
  y1 <- y[is]
  X1 <- X[is,]
  r2_boots[b] <- r2fxn( x1, y1, X1 )$r.squared
 }
 list( r2=r2, pgs_w0=coef, r2_boots=r2_boots, r2_sd=sd(r2_boots),r2_pval=lm_out$coefficients["x",4])
}

# bootrap R2 diff
r2_diff_boot <- function(x0,x1,y,B=1e4,X){
 r21 <- r2fxn(x1,y,X)$r.squared
 r20 <- r2fxn(x0,y,X)$r.squared
 delta_boots <- numeric()
 for( b in 1:B ){
  is <- sample( 1:length(x0), replace=T )
  x0b <- x0[is]
  x1b <- x1[is]
  yb <- y[is]
  Xb <- X[is,]
  delta_boots[b] <- (r2fxn(x0b,yb,Xb)$r.squared-r2fxn(x1b,yb,Xb)$r.squared) / r2fxn(x1b,yb,Xb)$r.squared
 }
 pv_oneside <- (1 + sum(delta_boots <= 0))/ (1+B)
 pval    <- 2*min(pv_oneside,1 - pv_oneside)
 list( r20=r20, r21=r21, ci_95=quantile(delta_boots,c(0.025,.975)), pval=pval)
}

r2_diff_log <- function(x0,x1,y,B=1e4,X){
 r21 <- r2fxn(x1,y,X)$r.squared
 r20 <- r2fxn(x0,y,X)$r.squared
 delta_boots <- numeric()
 for( b in 1:B ){
  is <- sample( 1:length(x0), replace=T )
  x0b <- x0[is]
  x1b <- x1[is]
  yb <- y[is]
  Xb <- X[is,]
  delta_boots[b] <- (r2fxn(x0b,yb,Xb)$r.squared-r2fxn(x1b,yb,Xb)$r.squared) / r2fxn(x0b,yb,Xb)$r.squared
 }
 pv_oneside <- (1 + sum(delta_boots <= 0))/ (1+B)
 pval    <- 2*min(pv_oneside,1 - pv_oneside)
 list( r20=r20, r21=r21, ci_95=quantile(delta_boots,c(0.025,.975)), pval=pval)
}

