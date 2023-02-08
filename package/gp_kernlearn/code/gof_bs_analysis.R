rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

save_slug <- 'package/gp_kernlearn/code/output/gof_bs/'
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
model.name <- 'model_Yrbf5-10_tr100_legd4'
saved.model.path <- paste0(save_slug, 'models/', model.name, '.RData')

# Either fit a model or load  a saved one.
if (is.na(model.name)) {
  deg <- 3
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

obs.mat <- matrix(unlist(model$Y.train), ncol=model$N.x, byrow = TRUE)
obs.cov.mat <- cov(obs.mat)
image(obs.cov.mat)
image(model$K.true)
image(model$K.est)

set.seed(2023)
norm.name <- 'frob_cov/'    # UPDATE 1/3
training.sumstat <- frob_norm_cov(obs.mat) # UPDATE 2/3
n.boot <- 2e2
boot.sumstats <- bootstrap_stat_dist(frob_norm_cov, model_data_gen, n.boot) # UPDATE 3/3

bs.output.path <- saved.model.path <- paste0(save_slug, norm.name, model.name)

png(paste0(bs.output.path, '_bs', n.boot,'dens.png'), height=6, width=6, units='in', res = 100)
plot(density(boot.sumstats), xlim=c(0, max(boot.sumstats, training.sumstat)))
points(training.sumstat, rep(0, length(training.sumstat)), col='red', pch=16)
dev.off()

bs.savepath <- paste0(bs.output.path, '_bs', n.boot,'.RData')
saveRDS(list(training.sumstat=training.sumstat,
             boot.sumstats=boot.sumstats),
        bs.savepath)

