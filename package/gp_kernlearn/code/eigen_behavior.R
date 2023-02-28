rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')
#source('package/gp_kernlearn/code/bootstrap.R')
source('package/gp_kernlearn/code/norms.R')

save_slug <- 'package/gp_kernlearn/code/output/gof_bs_exploration/'
dir.create(file.path(save_slug))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Training Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2023)

N.x <- 100         # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)
K <- rbf_kernel(X, 5, 0.05, 0.1)
kern.name <- 'rbf5-005-fz01'

true.sumstat <- eigen(K)$values
training.sumstat.list <- list()

n.train <- c(1e1, 1e2, 1e3, 1e4) 
for (i in 1:length(n.train)) {
  n.train.datasets <- n.train[i]
  
  Y.train <- list()
  for (r in 1:n.train.datasets) {
    err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), K))
    mu <- rep(0, N.x)
    Y <- mu + err
    Y.train[[r]] <- Y
  }
  obs.mat <- matrix(unlist(Y.train), ncol=N.x, byrow = TRUE)
  obs.cov.mat <- cov(obs.mat)
  training.sumstat.list[[i]] <- eigenvals_cov(obs.mat) # UPDATE 2/3
}

stat.dim <- 999
png(paste0(save_slug, 'n_datasets_eigenbehavior_2', kern.name,  '_nx', N.x,'.png'), height=6, width=6, units='in', res=100)
plot(log(ts(true.sumstat[1:stat.dim])), col='black',
     ylim=c(-6,6),
     lwd=2,
     xlab='Eigenvalue Order', ylab='log(Eigenvalue)', 
     main=NA)

colors = c('red', 'orange', 'green', 'blue', 'purple', 'magenta')
for (i in 1:length(n.train)) {
  lines(log(ts(training.sumstat.list[[i]][1:stat.dim])), col=colors[i], lty=2, lwd=2)  
}

legend('topright', legend=c(n.train, 'True'), col=c(colors[1:length(n.train)], 'black'), lty=c(rep(2, length(n.train)), 1), lwd=3)
dev.off()

saveRDS(list(training_sumstat_list=training.sumstat.list,
        true_sumstat=true.sumstat), paste0(save_slug, 'sumstats', kern.name,  '_nx', N.x,'.RData'))
