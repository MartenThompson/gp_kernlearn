rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Method Def ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
learn_kern <- function(Y, B.0, n.mcmc=1000) {
  B.0.svd <- svd(B.0)
  P.0 <- B.0.svd$u
  D.0 <- diag(B.0.svd$d)
  Q.0 <- B.0.svd$v
  #sum(B.0 - P.0%*%D.0%*%t(Q.0))
  
  H <- n.mcmc
  br.out <- stan_glm(Y ~ B.0-1, family=gaussian(), iter=H)
  beta.post <- matrix(NA, nrow=H, ncol=J)
  for (j in 1:J) {
    beta.post[,j] <- unlist(br.out$stanfit@sim$samples[[1]][[j]])
  } 
  
  beta.M2.post <- matrix(0, J, J)
  for (i in 1:H) {
    beta.M2.post <- beta.M2.post + t(t(beta.post[i,]))%*%beta.post[i,]
  }
  beta.M2.post <- beta.M2.post/H
  V.0 <- beta.M2.post
  
  V.1 <- D.0%*%t(Q.0)%*%V.0%*%Q.0%*%D.0
  V.1.eig <- eigen(V.1)
  P.1 <- V.1.eig$vectors
  D.1 <- diag(V.1.eig$values)
  
  # Final Estimators 
  Gamma.hat <- P.0%*%P.1
  Lambda.hat <- D.1
  Sigma.hat <- Gamma.hat%*%Lambda.hat%*%t(Gamma.hat)
  
  return(list(
    P.0=P.0, D.0=D.0, Q.0=Q.0,
    stan.fit=br.out, beta.post.samples=beta.post,
    V.0=beta.M2.post,
    V.1=V.1, P.1=P.1, D.1=D.1,
    Gamma.hat=Gamma.hat, Lambda.hat=Lambda.hat, Sigma.hat=Sigma.hat
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
X <- matrix(seq(-1,1, length.out=n), ncol=1)
#Sigma.true <- cubic_kernel(X, 1,1,1,1,1e-6)
Sigma.true <- rbf_kernel(X, 5, 0.05, 1)
Gamma.true <- eigen(Sigma.true)$vectors
Lambda.true <- diag(eigen(Sigma.true)$values)

set.seed(202)
Y <- t(mvtnorm::rmvnorm(1, sigma=Sigma.true))
Y <- Y-mean(Y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
J <- 9
basis_maker <- make_legendre1D_basis_maker(J-1)
#basis_maker <- make_monomial1D_basis_maker(J)
B.0 <- basis_maker(X)
plot(ts(B.0[,1]), ylim=c(-2,2))
for (j in 2:J) {
  lines(ts(B.0[,j]))
}

kl.out <- learn_kern(Y, B.0)

apply(kl.out$beta.post, 2, mean)
par(mfrow=c(2,2))
for (j in 1:J) {
  plot(ts(kl.out$beta.post[,j]))
}
par(mfrow=c(1,1))
#acf(beta.post)
plot(X, Y, pch=16)
points(X, kl.out$stan.fit$fitted.values, col='red', pch=16)
legend('topleft', legend=c('Observed', 'Reg Fit'), fill=c('black', 'red'))


# all the variance in direction of 1st eigenvector?
diag(kl.out$Lambda.hat)
diag(Lambda.true)[1:J]


image(Sigma.true)
image(kl.out$Sigma.hat)

par(mfrow=c(2,2))
for (j in 1:J) {
  plot(eigen(Sigma.true)$vectors[,j], kl.out$Gamma.hat[,j])  
}
par(mfrow=c(1,1))


plot(ts(Y), lwd=2)
lines(-ts(kl.out$Gamma.hat[,1]*sqrt(diag(kl.out$Lambda.hat)[1])), col='red', lwd=2)
lines(ts(kl.out$stan.fit$fitted.values), lty=2, lwd=2)
legend('topleft', legend=c('True', 'Reg Fit', 'Gamma 1'), col=c('black', 'black', 'red'), lty=c(1,2,1), lwd=2)

plot(kl.out$Sigma.hat[,1], kl.out$Sigma.hat[,2])
plot(Y,kl.out$Sigma.hat[,1])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Full Frequentist ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Not productive. Do not run.
# s2 <- 1
# beta.cov.freq <- s2*solve(t(B.0)%*%B.0)
# beta.hat <- solve(t(B.0)%*%B.0)%*%t(B.0)%*%Y
# V.0.freq <- beta.cov.freq + beta.hat%*%t(beta.hat)
# 
# V.1.freq <- D.0%*%t(Q.0)%*%V.0.freq%*%Q.0%*%D.0
# V.1.freq.eig <- eigen(V.1.freq)
# P.1.freq <- V.1.freq.eig$vectors
# D.1.freq <- diag(V.1.freq.eig$values)
# 
# 
# Gamma.hat.freq <- P.0%*%P.1.freq
# Lambda.hat.freq <- D.1.freq
# Sigma.hat.freq <- Gamma.hat.freq%*%Lambda.hat.freq%*%t(Gamma.hat.freq)
# 
# diag(D.1.freq)
# image(Sigma.hat.freq)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Replicated Datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
learn_kern_from_many <- function(B.0, V.0) {
  B.0.svd <- svd(B.0)
  P.0 <- B.0.svd$u
  D.0 <- diag(B.0.svd$d)
  Q.0 <- B.0.svd$v

  V.1 <- D.0%*%t(Q.0)%*%V.0%*%Q.0%*%D.0
  V.1.eig <- eigen(V.1)
  P.1 <- V.1.eig$vectors
  D.1 <- diag(V.1.eig$values)
  
  # Final Estimators 
  Gamma.hat <- P.0%*%P.1
  Lambda.hat <- D.1
  Sigma.hat <- Gamma.hat%*%Lambda.hat%*%t(Gamma.hat)
  
  return(list(
    P.0=P.0, D.0=D.0, Q.0=Q.0,
    V.1=V.1, P.1=P.1, D.1=D.1,
    Gamma.hat=Gamma.hat, Lambda.hat=Lambda.hat, Sigma.hat=Sigma.hat
  ))
}

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
R <- 200 # R = 200 has very good agreement b/t D.1 and Lambda
n <- 100
X <- matrix(seq(-1,1, length.out=n), ncol=1)
Sigma.true <- rbf_kernel(X, 5, 0.05, 1) # cubic_kernel(X, 1,1,1,1,1e-6)
Y.train <- gen_data(R, Sigma.true)

J <- 9
basis_maker <- make_legendre1D_basis_maker(J-1)
B.0 <- basis_maker(X)

old.model <- fit_kernlearn(B.0, basis_maker, Y.train) # R*10 seconds

# Results when V.0 = 2nd Moment(beta|Y) 
old.model$E.est # E[beta|Y]. not quite zero
beta.M2 <- old.model$V.est + t(t(old.model$E.est))%*%old.model$E.est
kl.out <- learn_kern_from_many(B.0, beta.M2)
image(kl.out$Sigma.hat)
range(kl.out$Sigma.hat)

diag(kl.out$Lambda.hat)
eigen(Sigma.true)[[1]][1:J]


# Nearly identical results when V.0 = Var(beta|Y) 
beta.cov <- old.model$V.est
kl.out <- learn_kern_from_many(B.0, beta.cov)
image(kl.out$Sigma.hat)
range(kl.out$Sigma.hat)

diag(kl.out$Lambda.hat)
eigen(Sigma.true)[[1]][1:J]




