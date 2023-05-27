rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/acv_est.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Gen ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(212)
R <- 20
n <- 50
X <- matrix(seq(-1,1, length.out=n), ncol=1)
sig.0 <- 1
sig.1 <- 1
sig.2 <- 1
fuzz <- 1e-5
#Sigma.true <- rbf_kernel(X, sig.0, sig.1, fuzz) 
#Sigma.true <- cubic_kernel(X, sig.0, sig.1,sig.2,sig.3,fuzz)
Sigma.true <- matern_kernel_euc(X, sig.0, sig.1, sig.2, fuzz)
Y.train <- gen_data(R, Sigma.true)

#plot(ts(Sigma.true[1,]))
#plot_yobs(X, Y.train, 5)

#data.name <- paste0('Ytrain_rbf',sig.0, '-', sig.1,'-',fuzz,'_R',R,'Nx',n,'/')
#data.name <- paste0('Ytrain_cubic',sig.0, '-', sig.1,'-', sig.2,'-', sig.3,'-',fuzz,'_R',R,'Nx',n,'/')
data.name <- paste0('Ytrain_matern',sig.0, '-', sig.1,'-', sig.2,'-',fuzz,'_R',R,'Nx',n,'/')
save.slug <- paste0('package/gp_kernlearn/code/output/acv/', data.name)
dir.create(save.slug)

saveRDS(Y.train, paste0(save.slug, 'Ytrain.RData'))

png(paste0(save.slug, paste0('Ytrain.png')), height=5, width=5, units='in', res = 100)
plot_yobs(X, Y.train, 3)
dev.off()

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


J <- 4
basis_maker <- make_legendre1D_basis_maker(J-1)
B.0 <- basis_maker(X)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis: Within & B/t Var ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#lag <- ifelse(n < 100 || R < 100, min(floor(0.75*R), floor(0.5*n)), 25)
lag <- min(floor(0.75*R), floor(0.5*n), 75)
acv.cov.within <- est_acv_combo(Y.train, B.0, 1e3, lag, 'covariance')
saveRDS(acv.cov.within, paste0(save.slug, 'acv_cov','_J', J, '_L', lag, '.RData'))
#acv.cor.within <- est_acv_within(Y.train, B.0, 1e3, lag, 'correlation')

# looks reasonable
#image(make_sigmahat_from_acv(acv.cov.within$acv.all[1,]))

# tightly correlated components, can probably use element-wise quantiles for HPD
#plot(acv.cov.within$acv.all[,2], acv.cov.within$acv.all[,3], pch=16)

alpha <- 0.05

png(paste0(save.slug, 'leg', 'J', J, '_L', lag, '_alpha', alpha, '_within.png'), height=5, width=5, units='in', res = 100)
plot_results(acv.cov.within$acv.all, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance', 
             '',
             Sigma.true[1,1], fuzz,
             alpha = alpha
             )
dev.off()


png(paste0(save.slug, 'leg', 'J', J, '_L', lag, '_alpha', alpha,'_within_between.png'), height=5, width=5, units='in', res = 100)
plot_results(acv.cov.within$acv.all.w.bt, 
             acf(unlist(Y.train), lag=lag, plot=FALSE, type='covariance')$acf, 
             Sigma.true[1,],
             'Covariance',
             '',
             Sigma.true[1,1], fuzz,
             alpha=alpha
             )
dev.off()
