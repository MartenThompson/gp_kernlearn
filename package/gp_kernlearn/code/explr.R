rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
X <- matrix(seq(-1,1, length.out=n), ncol=1)
#Sigma.true <- cubic_kernel(X, 1,1,1,1,1e-6)
Sigma.true <- rbf_kernel(X, 5, 0.05, 1e-6)
Gamma.true <- eigen(Sigma.true)$vectors
Lambda.true <- diag(eigen(Sigma.true)$values)

set.seed(202)
Y <- t(mvtnorm::rmvnorm(1, sigma=Sigma.true))
Y <- Y-mean(Y)
#plot(X,Y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Method ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

J <- 9
basis_maker <- make_legendre1D_basis_maker(J-1)
B.0 <- basis_maker(X)
plot(ts(B.0[,1]), ylim=c(-2,2))
for (j in 2:J) {
  lines(ts(B.0[,j]))
}
B.0.svd <- svd(B.0)
P.0 <- B.0.svd$u
D.0 <- diag(B.0.svd$d)
Q.0 <- B.0.svd$v
#sum(B.0 - P.0%*%D.0%*%t(Q.0))

H <- 1000
br.out <- stan_glm(Y ~ B.0-1, family=gaussian())
beta.post <- matrix(NA, nrow=H, ncol=J)
for (j in 1:J) {
  beta.post[,j] <- unlist(br.out$stanfit@sim$samples[[1]][[j]])
} 

apply(beta.post, 2, mean)
par(mfrow=c(2,2))
for (j in 1:J) {
  plot(ts(beta.post[,j]))
}
par(mfrow=c(1,1))
#acf(beta.post)
plot(X, Y, pch=16)
points(X, br.out$fitted.values, col='red', pch=16)

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Final Estimators ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gamma.hat <- P.0%*%P.1
Lambda.hat <- D.1
Sigma.hat <- Gamma.hat%*%Lambda.hat%*%t(Gamma.hat)


diag(Lambda.hat)
diag(Lambda.true)[1:J]

image(Sigma.true)
image(Sigma.hat)

par(mfrow=c(2,2))
for (j in 1:J) {
  plot(eigen(Sigma.true)$vectors[,j], Gamma.hat[,j])  
}
par(mfrow=c(1,1))



#V.temp <- B.0%*%beta.M2.post%*%t(B.0)
#V.temp.eig <- eigen(V.temp)
#zapsmall(V.temp.eig$values)
