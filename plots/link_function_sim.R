rm( list=ls() )

savefile <- 'sim.Rdata'
if( file.exists(savefile) ){
  load( savefile )
} else {

Ns    <- c(1e3,1e4,1e5)
S     <- 1e2
prev  <- .2
h2    <- .2
sig2E <- .2

nit <- 1e5
out <- array( NA, dim=c(3,4,2,nit) )
}

for( it in 1:nit )
  for( n in 1:3 )
{
  if( !is.na(out[n,1,1,it]) ) next

G   <- scale( sapply( runif(S,.05,.5), function(maf) rbinom(Ns[n],2,maf) ) )
E   <- scale( rbinom(Ns[n],1,.5) )

betas <- sqrt(h2/S) * rnorm(S)
PS  <- G %*% betas
z   <- PS + sqrt(sig2E) * E + sqrt(1-h2-sig2E)*rnorm(Ns[n])
y   <- as.numeric( z > quantile(z,1-prev) )


mains <- c( 'OLS on Liability', 'OLS on 0/1', 'Logit', 'Probit' )
out[n,1,1,it] <- summary( lm(   z ~ PS * E                                   ) )$coef['PS:E',4] 
out[n,2,1,it] <- summary( lm(   y ~ PS * E                                   ) )$coef['PS:E',4] 
out[n,3,1,it] <- summary( glm(  y ~ PS * E, family = binomial(link="logit")  ) )$coef['PS:E',4] 
out[n,4,1,it] <- summary( glm(  y ~ PS * E, family = binomial(link="probit") ) )$coef['PS:E',4] 

G1 <- G[,1]
out[n,1,2,it] <- summary( lm(   z ~ G1 * E                                   ) )$coef['G1:E',4] 
out[n,2,2,it] <- summary( lm(   y ~ G1 * E                                   ) )$coef['G1:E',4] 
out[n,3,2,it] <- summary( glm(  y ~ G1 * E, family = binomial(link="logit")  ) )$coef['G1:E',4] 
out[n,4,2,it] <- summary( glm(  y ~ G1 * E, family = binomial(link="probit") ) )$coef['G1:E',4] 


pdf( 'sim.pdf',h=6,w=18 )
par( mfcol=c(2,6), mar=c(5,5,1,1) )

breaks  <- seq(0-1e-8,1+1e-8,len=21)
mids    <- (breaks[-1]+breaks[-21])/2
xlabs   <- c( 'PRSxE p-value', 'SNPxE p-value' )
ylabs   <- c( 'PRSxE -log10(p)', 'SNPxE -log10(p)' )

for( k in 1:2 )
for( n in 1:3 ) 
try({

  plot( 0:1, c(0,20), main='', xlab=xlabs[k], ylab='Density', type='n' ) 
  for( j in 1:4 ){
    ys <- hist( out[n,j,k,], breaks=breaks, plot=F )$density
    lines( mids, ys, col=j, lwd=2 )
  }
  abline( h=1, col='grey', lwd=2, lty=2 ) 
  legend( 'topright', fill=c(1:4,'grey'), leg=c(mains,'Null'), bty='n' )
  legend( 'right', leg=paste0( 'N=', Ns[n] ), bty='n' )

  lims <- c(0,c(10,6)[k])
  plot( lims, lims, type='n', ylab=ylabs[k], xlab='Expected -log10(p)' )

  my_qq <- function(y,...){
    n <- sum(!is.na(y))
    y <- c(0,sort(-log10(y)))
    x <- c(0,sort(-log10(1:n/(1+n))))
    points( x, y, pch=16,... )
    lines ( x, y, lwd=2 ,... )
  }
  for( j in 1:4 )
    my_qq( out[n,j,k,], col=j )
  abline(a=0,b=1)
  legend( 'right', leg=paste0( 'N=', Ns[n] ), bty='n' )

})
dev.off()

save.image( file=savefile )

}
