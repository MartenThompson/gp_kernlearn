rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/acv_est.R')

gen_data <- function(n.train.datasets, Kern) {
  N.x <- ncol(Kern)
  Y.train <- list()
  for (r in 1:n.train.datasets) {
    err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), Kern))
    mu <- rep(0, N.x)
    Y <- mu + err
    Y.train[[r]] <- Y
  }  
  return(Y.train)
}


set.seed(212)
R <- 100 # R = 200 has very good agreement b/t D.1 and Lambda
n <- 100
X <- matrix(seq(-1,1, length.out=n), ncol=1)
Sigma.true <- rbf_kernel(X, 5, 0.05, 0.1) # cubic_kernel(X, 1,1,1,1,1e-6)
Y.train <- gen_data(R, Sigma.true)

J <- 7
basis_maker <- make_legendre1D_basis_maker(J-1)
#basis_maker <- make_monomial1D_basis_maker(J)
B.0 <- basis_maker(X)

acv.out <- est_sigma_acv(Y.train, B.0, 1e3)
dim(acv.out$acv.all)
#image(acv.out$acv.all)
image(make_sigmahat_from_acv(acv.out$acv.all[3000,]))


plot(acv.out$acv.all[,2],acv.out$acv.all[,3], pch=16)

plot(ts(acv.out$acv.all[1,]), col=rgb(0,0,0,0.15), ylim=c(-5,15))
for (i in 2:nrow(acv.out$acv.all)) {
  if (i %% 500 == 0) {
    lines(ts(acv.out$acv.all[i,]), col=rgb(0,0,0,0.15))  
  }
}
alpha <- 0.05
pd.quant.pointwise <- apply(acv.out$acv.all, 2, function(c) {quantile(c, c(alpha/2, .5, 1-alpha/2))})
polygon(c(1:100, 100:1), c(pd.quant.pointwise[1,], rev(pd.quant.pointwise[3,])), border=NA, col=rgb(1,0,0,0.25))
lines(ts(pd.quant.pointwise[2,]), col='red', lwd=2)
lines(ts(Sigma.true[1,]), col='green', lwd=2)


Y.mat <- matrix(unlist(Y.train), ncol=n, byrow = TRUE)
Y.acv <- acf(Y.mat, lag.max = 99, type='covariance', plot=FALSE)
plot(ts(Y.acv$acf[,1,1]))
for (i in 2:100) {
  lines(ts(Y.acv$acf[,i,i]))
}

