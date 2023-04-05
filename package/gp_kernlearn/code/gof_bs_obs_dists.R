rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/norms.R')
source('package/gp_kernlearn/code/kernlearn.R')

save_slug <- 'package/gp_kernlearn/code/output/gof_bs/obs_dists/'
dir.create(file.path(save_slug))


gen_data <- function(n.train.datasets, N.x, Kern) {
  Y.train <- list()
  for (r in 1:n.train.datasets) {
    err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), Kern))
    mu <- rep(0, N.x)
    Y <- mu + err
    Y.train[[r]] <- Y
  }
  obs.mat <- matrix(unlist(Y.train), ncol=N.x, byrow = TRUE)
  return(obs.mat)  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Training Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2023)

replications <- 1e2
N.x <- 100         # number of x locations
n.train.datasets <- 100 
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)
 
#L <- 5
#norm.name <- paste0('wem_neg2pow_L', L) # UPDATE 1/2
#norm_func <- weighted_eval_cov_maker(L=L, schedule = rep(2,L)^c(-L:-1))# UPDATE 2/2
norm.name <- 'whichmax_eigencentraldiff'
norm_func <- max_centraldifflogeigen_cov

obs.sumstats <- data.frame(
  lin=rep(NA, replications),
  lin_nofz=rep(NA, replications),
  quad=rep(NA, replications),
  quad_nofz=rep(NA, replications),
  cubic=rep(NA, replications),
  cubic_nofz=rep(NA, replications),
  rbf5_05=rep(NA, replications),
  rbf5_05_nofz=rep(NA, replications),
  rbf5_005=rep(NA, replications),
  rbf5_005_nofz=rep(NA, replications)
)

for (i in 1:replications) {
  cat(i, ' ')
  K <- lin_kernel(X, 1, 1, 0.1)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$lin[i] <- norm_func(obs.mat)
  
  K <- lin_kernel(X, 1, 1, 0)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$lin_nofz[i] <- norm_func(obs.mat)
  
  K <- quad_kernel(X, 1, 1, 1, 0.1)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$quad[i] <- norm_func(obs.mat)
  
  K <- quad_kernel(X, 1, 1, 1, 0)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$quad_nofz[i] <- norm_func(obs.mat)
  
  K <- cubic_kernel(X, 1, 1, 1, 1, 0.1)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$cubic[i] <- norm_func(obs.mat)
  
  K <- cubic_kernel(X, 1, 1, 1, 1, 0)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$cubic_nofz[i] <- norm_func(obs.mat)
  
  K <- rbf_kernel(X, 5, 0.5, 0.1)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$rbf5_05[i] <- norm_func(obs.mat)
  
  K <- rbf_kernel(X, 5, 0.5, 0)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$rbf5_05_nofz[i] <- norm_func(obs.mat)
  
  K <- rbf_kernel(X, 5, 0.05, 0.1)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$rbf5_005[i] <- norm_func(obs.mat)
  
  K <- rbf_kernel(X, 5, 0.05, 0)
  obs.mat <- gen_data(n.train.datasets, N.x, K)
  obs.sumstats$rbf5_005_nofz[i] <- norm_func(obs.mat)
}
summary(obs.sumstats)

write.csv(obs.sumstats, paste0(save_slug, norm.name, '_sumstats.csv'), row.names = FALSE)

#png(paste0(save_slug, norm.name,'_dens.png'), height=6, width=6, units='in', res = 100)
plot(density(obs.sumstats$lin), lwd=2, col='red',
     xlim=c(min(obs.sumstats), max(obs.sumstats)),
     ylim=c(0, max(density(obs.sumstats$rbf5_005)$y +0.005, density(obs.sumstats$lin)$y +0.005)),
     xlab=NA, main=NA)
lines(density(obs.sumstats$lin_nofz), lwd=2, col='red', lty=2)
lines(density(obs.sumstats$quad), lwd=2, col='orange')
lines(density(obs.sumstats$quad_nofz), lwd=2, col='orange', lty=2)
lines(density(obs.sumstats$cubic), lwd=2, col='green')
lines(density(obs.sumstats$cubic_nofz), lwd=2, col='green', lty=2)
lines(density(obs.sumstats$rbf5_05), lwd=2, col='blue')
lines(density(obs.sumstats$rbf5_05_nofz), lwd=2, col='blue', lty=2)
lines(density(obs.sumstats$rbf5_005), lwd=2, col='magenta')
lines(density(obs.sumstats$rbf5_005_nofz), lwd=2, col='magenta', lty=2)
legend('topright', 
       c('Linear', 'Quadratic', 'Cubic', 'RBF Smooth', 'RBF Rough'), 
       col=c('red', 'orange', 'green', 'blue', 'magenta'), lty=rep(1, 5), lwd=3)
#dev.off()

require(ggplot2)
require(reshape2)

plot.data <- melt(obs.sumstats[,c(1,3,5,7)])

ggplot(data = plot.data) +
  geom_histogram(aes(x=value, group=variable, fill=variable), position = 'dodge', bins=5)


plot.data <- melt(obs.sumstats[,c(7,9)])
ggplot(data=plot.data) + 
  geom_boxplot(aes(x=value, group=variable))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### RBF Exploration ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rbf.vars <- c(1,5,10)
rbf.ls <- c(5, 0.5, 0.05, 0.005)
rbf.fz <- 0.1
rbf.params <- expand.grid(rbf.vars, rbf.ls, KEEP.OUT.ATTRS = F)
datagen.names <- apply(rbf.params, 1, function(x) {
  paste0('rbf',x[1], '-', x[2], '-fz', rbf.fz)
})

set.seed(2023)

replications <- 1e2
N.x <- 100         # number of x locations
n.train.datasets <- 100 
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)

norm.name <- 'whichmax_eigencentraldiff_rbf_fz'
norm_func <- max_centraldifflogeigen_cov

obs.sumstats <- data.frame(matrix(NA, nrow=replications, ncol=length(datagen.names)))
colnames(obs.sumstats) <- datagen.names


for (i in 1:replications) {
  cat(i, ' ')
  for (j in 1:nrow(rbf.params)) {
    v <- rbf.params[j,1]
    l <- rbf.params[j,2]
    K <- rbf_kernel(X, v, l, rbf.fz)
    obs.mat <- gen_data(n.train.datasets, N.x, K)
    obs.sumstats[i,j] <- norm_func(obs.mat)
  }
}

write.csv(obs.sumstats, paste0(save_slug, norm.name, '_sumstats.csv'), row.names = FALSE)

plot.data <- melt(obs.sumstats)


dg.labels <- lapply(strsplit(levels(plot.data$variable), '-'), function(l) {
  paste0('ls=', l[2], ', var=', substr(l[1], 4, nchar(l[1])))})


png(paste0(save_slug, norm.name,'_dens.png'), height=6, width=6, units='in', res = 100)
ggplot(plot.data) + 
  geom_boxplot(aes(x=variable, y=value)) +
  scale_x_discrete(name='Data Generating Process', labels=dg.labels) + 
  scale_y_continuous(name='Index of Max Central Difference') + 
  coord_flip() +
  theme_bw()
dev.off()

#K <- rbf_kernel(X, 5, 5, 0)
#obs.mat <- gen_data(n.train.datasets, N.x, K)
#lot(obs.mat[1,])


K <- matern_kernel_euc(X, 1, 10)
obs.mat <- gen_data(n.train.datasets, N.x, K)
plot(ts(log(eigen(cov(obs.mat))[[1]])))
plot(ts(central_difference(log(eigen(cov(obs.mat))[[1]]))))
  
plot(ts(eigen(cov(obs.mat))[[1]]))
plot(ts(central_difference(eigen(cov(obs.mat))[[1]])))

