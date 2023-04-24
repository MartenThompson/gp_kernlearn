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

plot_results <- function(acv.est.all, acv.ytrain, acv.true, ylab='cov', title='', alpha=0.05, thin.mod=500) {
  lag <- ncol(acv.est.all)
  ymin <- min(quantile(acv.est.all, 0.01), acv.ytrain, acv.true)
  ymax <- max(quantile(acv.est.all,0.99), acv.ytrain, acv.true)
  
  plot(ts(acv.est.all[1,]), col=rgb(0,0,0,0.15), 
       ylim=c(ymin, ymax), xaxs='i',
       xlab='lag', ylab=ylab, main=title)
  for (i in 2:nrow(acv.est.all)) {
    if (i %% thin.mod == 0) {
      lines(ts(acv.est.all[i,]), col=rgb(0,0,0,0.15))  
    }
  }
  
  pd.quant.pointwise <- apply(acv.est.all, 2, function(c) {quantile(c, c(alpha/2, .5, 1-alpha/2))})
  polygon(c(1:lag, lag:1), c(pd.quant.pointwise[1,], rev(pd.quant.pointwise[3,])), border=NA, col=rgb(1,0,0,0.25))
  lines(ts(pd.quant.pointwise[2,]), col='red', lwd=2)
  lines(ts(acv.true), col='green', lwd=2)
  lines(ts(acv.ytrain), col='blue', lwd=2, lty=2)
  legend('topright', legend=c('True', 'Median Est', 'Y.train'), col=c('green', 'red', 'blue'), lwd=3, lty=c(1,1,2))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Gen ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(212)
R <- 50
n <- 100
X <- matrix(seq(-1,1, length.out=n), ncol=1)
sig.0 <- 1
sig.1 <- 0.1
fuzz <- 1e-1
Sigma.true <- rbf_kernel(X, sig.0, sig.1, fuzz) 
Y.train <- gen_data(R, Sigma.true)

plot(Y.train[[1]], col=rgb(1,1,1), ylim=c(min(unlist(Y.train)), max(unlist(Y.train))))
for (i in 1:5) {
  lines(Y.train[[i]])
}

# Causes of underestimating variance:
  # Noise: think of a flat line fit (zero var), and obs noisey (some var).
  # x window
  # var within and var b/t
  # J small

# both within and b/t variance
var(unlist(Y.train))
# just within
temp <- lapply(Y.train, function(x) {var(x)})
summary(unlist(temp))


J <- 8
basis_maker <- make_legendre1D_basis_maker(J-1)
B.0 <- basis_maker(X)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis: Within Funct Var Only ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lag <- min(n-1, 25)
acv.cov.within <- est_acv_within(Y.train, B.0, 1e3, lag, 'covariance')
#acv.cor.within <- est_acv_within(Y.train, B.0, 1e3, lag, 'correlation')

# looks reasonable
image(make_sigmahat_from_acv(acv.cov.within$acv.all[1,]))

# tightly correlated components, can probably use element-wise quantiles for HPD
#plot(acv.cov.within$acv.all[,2], acv.cov.within$acv.all[,3], pch=16)

plot_results(acv.cov.within$acv.all, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). Within. J=', J, '. Nx=', n, '. R=', R),
             )


plot_results(acv.cov.within$acv.all.w.bt, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). Within. J=', J, '. Nx=', n, '. R=', R),
)



plot(ts(acf(Y.train[[1]], lag=15, plot=FALSE, type='covariance')$acf), ylim=c(-4,8))
for (i in 2:length(Y.train)) {
  lines(ts(acf(Y.train[[i]], lag=15, plot=FALSE, type='covariance')$acf))
}




plot_results(temp, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). Within. J=', J),
)


plot_results(acv.cor.within$acv.all, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='correlation')$acf, 
             cov2cor(Sigma.true)[1,],
             'Correlation', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). Within. J=', J))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis: Within & Between Funct Var ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
acv.cov.bt <- est_acv_bt(Y.train, B.0, 1e3, lag, 'covariance')
acv.cor.bt <- est_acv_bt(Y.train, B.0, 1e3, lag, 'correlation')

plot_results(acv.cov.bt$acv.all, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). B/t. J=', J))

plot_results(acv.cor.bt$acv.all, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='correlation')$acf, 
             cov2cor(Sigma.true)[1,],
             'Correlation', paste0('RBF(', sig.0,',', sig.1, ',', fuzz, '). B/t. J=', J))

###


temp.stan <- stan_beta_post(Y.train[[1]], B.0, 1e3)

p <- ncol(B.0)
df <- n-p
V.beta <- solve(t(B.0)%*%B.0)
beta.ols <- solve(t(B.0)%*%B.0)%*%t(B.0)%*%Y.train[[1]]
s2 <- ((1/(n-p))*t(Y.train[[1]] - B.0%*%beta.ols)%*%(Y.train[[1]] - B.0%*%beta.ols))[1,1]
t.scale.mat <- s2*(matrix(1, ncol=n, nrow=n) + B.0%*%V.beta%*%t(B.0))
t.pp.var <- (df/(df-2))*t.scale.mat


