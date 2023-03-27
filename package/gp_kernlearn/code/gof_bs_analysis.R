rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/stationary.R')

save_slug <- 'package/gp_kernlearn/code/output/gof_bs/'
dir.create(file.path(save_slug))
dir.create(file.path(paste0(save_slug, 'models')))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Training Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2023)
N.x <- 100         # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)
#K <- lin_kernel(X, 1, 1, 0.1)
#K <- quad_kernel(X, 1, 1, 1, 0.1)
#K <- cubic_kernel(X, 1, 1, 1, 1, 0.1)
K <- rbf_kernel(X, 5, 0.05, 0.1) 
gen.name <- 'rbf5-005-fz01' # 'lin1-1' 
n.train.datasets <- 100 # this actually needs to be 1e3 before cov(Y.train) looks like K

Y.train <- list()
for (r in 1:n.train.datasets) {
  err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), K))
  mu <- rep(0, N.x)
  Y <- mu + err
  Y.train[[r]] <- Y
}

plot(Y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Fit/Load Model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model.name <- 'model_Ycubic_tr100_legd3' # NA
saved.model.path <- paste0(save_slug, 'models/', model.name, '.RData')

# Either fit a model or load  a saved one.
if (is.na(model.name)) {
  deg <- 2
  leg_basis_maker <- make_legendre1D_basis_maker(degree = deg)
  b.X <- leg_basis_maker(X)
  model <- fit_kernlearn(b.X, leg_basis_maker, Y.train)
  model$K.true <- K
  model$N.x <- N.x
  model$n.train.datasets <- n.train.datasets
  saveRDS(model, paste0(save_slug, 'models/model_', 'Y', gen.name, 
                        '_tr', n.train.datasets, '_legd', deg ,'.RData'))
} else {
  model <- readRDS(saved.model.path)  
  
  model_data_gen_single <- function() {
    K.est <- model$K.est
    return(matrix(mvtnorm::rmvnorm(1,sigma=K.est), ncol=1))
  }
  X <- model$b.X.train[,2]
  s.out.mon3 <- make_stationary(model$K.est, X, dist_euc, modfit_monopoly_maker(3))
  model_data_gen_stationary_single <- function() {
    return(matrix(mvtnorm::rmvnorm(1,sigma=s.out.mon3$K.stationary), ncol=1))
  }
}


plot(model_data_gen_single())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Bootstrapping ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('package/gp_kernlearn/code/bootstrap.R')
source('package/gp_kernlearn/code/norms.R')

model_data_gen <- function() {
  # matrix with n.train.datasets rows, each a synthetic obs from our model
  t(sapply(1:model$n.train.datasets, function(x){model_data_gen_single()}))
}
#model_data_gen_stationary <- function() {
  # matrix with n.train.datasets rows, each a synthetic obs from our model
#  t(sapply(1:model$n.train.datasets, function(x){model_data_gen_stationary_single()}))
#}

obs.mat <- matrix(unlist(model$Y.train), ncol=model$N.x, byrow = TRUE)
#obs.cov.mat <- cov(obs.mat)
#image(obs.cov.mat)
#image(model$K.true)
#image(model$K.est)
#image(model$V.est)

set.seed(2023)
norm.name <- 'wem_uniform/'    # UPDATE 1/3
K = 5
norm <- weighted_eval_cov_maker(K=K, schedule = rep(1,K))
training.sumstat <- norm(obs.mat) # UPDATE 2/3
n.boot <- 2e2
boot.sumstats <- bootstrap_stat_dist(norm, model_data_gen, n.boot) # UPDATE 3/3
bs.output.path <- paste0(save_slug, norm.name, model.name, '_K', K)

png(paste0(bs.output.path, '_bs', n.boot,'dens.png'), height=4, width=4, units='in', res = 100)
stat.dim <- 1
plot(density(boot.sumstats[,stat.dim]), xlim=c(0, max(boot.sumstats[,stat.dim], training.sumstat[stat.dim])),
     xlab=NA, main=NA)
points(training.sumstat[stat.dim], c(0), col='red', pch=16)
dev.off()

bs.savepath <- paste0(bs.output.path, '_bs', n.boot,'.RData')
saveRDS(list(training.sumstat=training.sumstat,
             K = K, norm=norm,
             boot.sumstats=boot.sumstats),
        bs.savepath)


#dat.in <- readRDS('package/gp_kernlearn/code/output/gof_bs/eigenvals_cov/model_Ycubic_tr100_legd1_bs200.RData')
#training.sumstat <- dat.in$training.sumstat
#boot.sumstats <- dat.in$boot.sumstats

# stat.dim <- 40
# png(paste0(bs.output.path, '_bs', n.boot, '_eigenmatch', stat.dim, '.png'), height=6, width=6, units='in', res=100)
# plot(log(ts(eigen(model$K.true)$values[1:stat.dim])), col='green',
#      ylim=log(c(min(training.sumstat[1:stat.dim], apply(boot.sumstats[,1:stat.dim], 2, mean)),
#             max(training.sumstat[1:stat.dim], apply(boot.sumstats[,1:stat.dim], 2, mean)))),
#      lwd=2,
#      xlab='Eigenvalue Order', ylab='log(Eigenvalue)', 
#      main=model.name)
# lines(log(ts(training.sumstat[1:stat.dim])), col='blue', lwd=2)
# boot.means <- apply(boot.sumstats[,1:stat.dim], 2, mean)
# boot.sds <- apply(boot.sumstats[,1:stat.dim], 2, sd)
# lines(log(ts(boot.means)), col='magenta', lwd=2)
# polygon(c(1:stat.dim, rev(1:stat.dim)), c(log(boot.means + 2*boot.sds), rev(log(boot.means-2*boot.sds))),
#         col=rgb(1,0,1,0.25), border=NA)
# legend('topright', c('True', 'Obs', 'Mean Boot'), col=c('green', 'blue', 'magenta'), lty=1, lwd=3)
# dev.off()

# TODO: bootstrap inconsistent for model selection. Parametric bs is (and is it equiv to AIC, BIC depending on resample size). 
# Also, dependent data.
