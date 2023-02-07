rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/gof_bootstrap.R')


save_slug <- 'package/gp_kernlearn/code/output/gof_bs_frob/'
dir.create(file.path(save_slug))
dir.create(file.path(paste0(save_slug, '/models')))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Training Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2023)
N.x <- 100         # number of x locations
X <- matrix(seq(-5,5,length.out=N.x), nrow=N.x, ncol=1)
#K <- lin_kernel(X, 1, 1/25, 0.1)
#K <- quad_kernel(X, 1, 1/25, 1/25, 0.1)
#K <- cubic_kernel(X, 1, 1/25, 1/25, 1/25, 0.1)
K <- rbf_kernel(X, 5, 10)
gen.name <- 'rbf5-10'
n.train.datasets <- 20

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
saved.model.path <- paste0(save_slug, 'models/model_Ycubic_tr20_legd3.RData')

# Either fit a model or load  a saved one.
if (is.na(saved.model.path)) {
  deg <- 3
  leg_basis_maker <- make_legendre1D_basis_maker(degree = deg)
  b.X <- leg_basis_maker(X)
  model <- fit_kernlearn(b.X, leg_basis_maker, Y.train)
  saveRDS(model, paste0(save_slug, 'models/model_', 'Y', gen.name, '_tr', n.train.datasets, '_legd', deg ,'.RData'))
} else {
  model <- readRDS(saved.model.path)  
  
  model_data_gen <- function() {
    K.est <- model$K.est
    return(matrix(mvtnorm::rmvnorm(1,sigma=K.est), ncol=1))
  }
}


plot(model_data_gen())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Bootstrapping ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frob.norm <- function(D) {
  return(norm(D, type='F'))
}

training.sumstats <- unlist(lapply(Y.train, frob.norm))
n.boot <- 1e3
boot.sumstats <- bootstrap_stat_dist(frob.norm, model_data_gen, n.boot)


plot(density(boot.sumstats), xlim=c(0, max(boot.sumstats, training.sumstats)))
points(training.sumstats, rep(0, length(training.sumstats)), col='red', pch=16)

bs.savepath <- paste0(strsplit(saved.model.path, '.', fixed=T)[[1]][1], '_bs.RData')
saveRDS(list(training.sumstats=training.sumstats,
             boot.sumstats=boot.sumstats),
        bs.savepath)

