setwd('~/Git/gp_kernlearn/')
library(rstanarm)
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

rbf_kernel <- function(X, alpha.0, alpha.1, phi=1e-6) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  if (p >= N) {
    stop("Have not thought through p > N.")
  }
  
  K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
  diag(K) <- diag(K) + phi
  return(K)
}

frobenius_norm <- function(mat) {
  sum(mat^2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### p = 1. Y= 0 + GP(0, rbf) rough ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
save_slug <- 'package/gp_kernlearn/vignettes/output/br_0mean_rbf_p1/'
dir.create(file.path(save_slug))

N <- 50
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- rbf_kernel(X, 1, 5, 0.1)
saveRDS(X, paste0(save_slug, 'X.RData'))
saveRDS(K, paste0(save_slug, 'K.RData'))
image(K, main='K True')

leg.deg <- 6
leg_basis_maker <- make_legendre1D_basis_maker(degree=leg.deg)
b.X <- leg_basis_maker(X)
dir.create(file.path(paste0(save_slug, 'leg', leg.deg, '/')))
saveRDS(b.X, paste0(save_slug, 'leg', leg.deg, '/bX.RData'))

reps <- 50
beta.samples.many <- NA
Y.history <- list()
err.history <- list()
set.seed(2022)
for (r in 1:reps) {
  cat('Replication', r, '\n')
  err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
  mu <- 0*(-2*X + 0.5*X^2 + 0.1*X^3) 
  Y <- mu + err
  Y.history[[r]] <- Y
  err.history[[r]] <- err
  
  #plot(X, Y, main=paste0('Rep ', r))
  #points(X, mu, col='green')
  #points(X, err, col='red')
  
  br.output <- stan_glm(Y ~ b.X[,2:(leg.deg+1)], family=gaussian())

  beta.samples <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[1]]), ncol=1)
  for (i in 2:(leg.deg+1)) {
    new <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[i]]), ncol=1)
    beta.samples <- cbind(beta.samples, new)
  }

  if (1==r) {
    beta.samples.many <- beta.samples
  } else {
    beta.samples.many <- rbind(beta.samples.many, beta.samples)
  }
    
}

saveRDS(Y.history, paste0(save_slug, 'leg', leg.deg, '/Yhist.RData'))
saveRDS(err.history, paste0(save_slug, 'leg', leg.deg, '/errhist.RData'))
saveRDS(beta.samples.many, paste0(save_slug, 'leg', leg.deg, '/betasamples.RData'))

E <- apply(beta.samples.many, 2, mean)
V <- cov(beta.samples.many)
saveRDS(E, paste0(save_slug, 'leg', leg.deg, '/E.RData'))
saveRDS(V, paste0(save_slug, 'leg', leg.deg, '/V.RData'))

#X.test <- matrix(sort(runif(100, -5, 5)), ncol=1)
post.ests <- posterior_test_meanvar_brev(E, V, X, Y.history[[3]], X.test, leg_basis_maker)

png(paste0(save_slug, 'leg', leg.deg, '/response_eg.png'))
plot(X, Y.history[[3]], ylim=c(-2,3))
#points(X, mu, col='green')
points(X.test, post.ests$post.test.mean, col='blue', pch=16)
points(rep(X.test,2), c(post.ests$post.test.mean + sqrt(diag(post.ests$kern.pieces$K.star)),
                        post.ests$post.test.mean - sqrt(diag(post.ests$kern.pieces$K.star))),
       pch=16, col=rgb(0,0,1,0.5))
dev.off()

png(paste0(save_slug, 'leg', leg.deg, '/K_train.png'))
image(post.ests$kern.pieces$K.train, main=paste0('Legendre ', leg.deg))
dev.off()
#image(post.ests$kern.pieces$K.star)

saveRDS(post.ests$kern.pieces$K.train, paste0(save_slug, 'leg', leg.deg, '/Ktrain_est.RData'))

#K[1:6,1:6]
#post.ests$kern.pieces$K.train[1:6,1:6]
#image(K)
#image(post.ests$kern.pieces$K.train + diag(0.5, nrow=dim(K)[1], ncol=dim(K)[1]))

norm <- frobenius_norm(K - post.ests$kern.pieces$K.train)
saveRDS(norm, paste0(save_slug, 'leg', leg.deg, '/norms.RData'))

#####
#####


analysis.dirs <- list.files(save_slug, pattern = 'leg*')
norms <- rep(NA, length(analysis.dirs))
for (i in 1:length(analysis.dirs)) {
  d <- analysis.dirs[i]
  norms[i] <- readRDS(paste0(save_slug, d, '/norms.RData'))
}
plot(NA, NA, xlim=c(1,10), ylim=c(0,200), xlab='Legendre Degree', ylab='Frobenius |K - Kest|')
lines(1:length(analysis.dirs), norms, lwd=2)




plot(NA,NA, xlim=c(-5,5), ylim=c(-3,3), xlab='X', ylab='Y')
lines(X, Y.history[[1]], lwd=2)
lines(X, Y.history[[2]], lwd=2, col='red')
lines(X, Y.history[[10]], lwd=2, col='blue')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### p = 1. Y=0 + GP(0, rbf) smooth ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
save_slug <- 'package/gp_kernlearn/vignettes/output/br_0mean_rbf_p1_2/'
dir.create(file.path(save_slug))

N <- 50
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- rbf_kernel(X, 1, 5, 0)
saveRDS(X, paste0(save_slug, 'X.RData'))
saveRDS(K, paste0(save_slug, 'K.RData'))

png(paste0(save_slug, 'K_true.png'))
image(K, main='K True')
dev.off()

leg.deg <- 10
leg_basis_maker <- make_legendre1D_basis_maker(degree=leg.deg)
b.X <- leg_basis_maker(X)
dir.create(file.path(paste0(save_slug, 'leg', leg.deg, '/')))
saveRDS(b.X, paste0(save_slug, 'leg', leg.deg, '/bX.RData'))

reps <- 50
beta.samples.many <- NA
Y.history <- list()
err.history <- list()
set.seed(2022)
for (r in 1:reps) {
  cat('Replication', r, '\n')
  err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
  mu <- 0*(-2*X + 0.5*X^2 + 0.1*X^3) 
  Y <- mu + err
  Y.history[[r]] <- Y
  err.history[[r]] <- err
  
  #plot(X, Y, main=paste0('Rep ', r))
  #points(X, mu, col='green')
  #points(X, err, col='red')
  
  br.output <- stan_glm(Y ~ b.X[,2:(leg.deg+1)], family=gaussian())
  
  beta.samples <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[1]]), ncol=1)
  for (i in 2:(leg.deg+1)) {
    new <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[i]]), ncol=1)
    beta.samples <- cbind(beta.samples, new)
  }
  
  if (1==r) {
    beta.samples.many <- beta.samples
  } else {
    beta.samples.many <- rbind(beta.samples.many, beta.samples)
  }
  
}

saveRDS(Y.history, paste0(save_slug, 'leg', leg.deg, '/Yhist.RData'))
saveRDS(err.history, paste0(save_slug, 'leg', leg.deg, '/errhist.RData'))
saveRDS(beta.samples.many, paste0(save_slug, 'leg', leg.deg, '/betasamples.RData'))

E <- apply(beta.samples.many, 2, mean)
V <- cov(beta.samples.many)
saveRDS(E, paste0(save_slug, 'leg', leg.deg, '/E.RData'))
saveRDS(V, paste0(save_slug, 'leg', leg.deg, '/V.RData'))

X.test <- matrix(sort(runif(100, -5, 5)), ncol=1)
post.ests <- posterior_test_meanvar_brev(E, V, X, Y.history[[3]], X.test, leg_basis_maker)

png(paste0(save_slug, 'leg', leg.deg, '/response_eg_leg', leg.deg, '.png'))
plot(X, Y.history[[3]], ylim=c(-2,3))
#points(X, mu, col='green')
points(X.test, post.ests$post.test.mean, col='blue', pch=16)
points(rep(X.test,2), c(post.ests$post.test.mean + sqrt(diag(post.ests$kern.pieces$K.star)),
                        post.ests$post.test.mean - sqrt(diag(post.ests$kern.pieces$K.star))),
       pch=16, col=rgb(0,0,1,0.5))
dev.off()

png(paste0(save_slug, 'leg', leg.deg, '/K_train_leg', leg.deg, '.png'))
image(post.ests$kern.pieces$K.train, main=paste0('Legendre ', leg.deg))
dev.off()
#image(post.ests$kern.pieces$K.star)

saveRDS(post.ests$kern.pieces$K.train, paste0(save_slug, 'leg', leg.deg, '/Ktrain_est.RData'))

#K[1:6,1:6]
#post.ests$kern.pieces$K.train[1:6,1:6]
#image(K)
#image(post.ests$kern.pieces$K.train + diag(0.5, nrow=dim(K)[1], ncol=dim(K)[1]))


norm <- frobenius_norm(K - post.ests$kern.pieces$K.train)
saveRDS(norm, paste0(save_slug, 'leg', leg.deg, '/norms.RData'))

#####
#####


analysis.dirs <- list.files(save_slug, pattern = 'leg*')
norms <- rep(NA, length(analysis.dirs))
for (i in 1:length(analysis.dirs)) {
  d <- analysis.dirs[i]
  norms[i] <- readRDS(paste0(save_slug, d, '/norms.RData'))
}
png(paste0(save_slug, 'frob_norm.png'))
plot(NA, NA, xlim=c(1,10), ylim=c(0,200), xlab='Legendre Degree', ylab='Frobenius |K - Kest|')
lines(1:length(analysis.dirs), norms, lwd=2)
dev.off()


png(paste0(save_slug, 'Ytrue_subset.png'))
plot(NA,NA, xlim=c(-5,5), ylim=c(-3,3), xlab='X', ylab='Y')
lines(X, Y.history[[1]], lwd=2)
lines(X, Y.history[[2]], lwd=2, col='red')
lines(X, Y.history[[10]], lwd=2, col='blue')
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### p = 1. Y=poly(x) + GP(0, rbf) smooth ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#mu <- 0*(-2*X + 0.5*X^2 + 0.1*X^3) 