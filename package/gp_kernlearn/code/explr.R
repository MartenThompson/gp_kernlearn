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
Sigma.true <- rbf_kernel(X, 5, 0.05, 1)
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
#basis_maker <- make_monomial1D_basis_maker(J)
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
br.out <- stan_glm(Y ~ B.0-1, family=gaussian(), iter=H)
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
legend('topleft', legend=c('Observed', 'Reg Fit'), fill=c('black', 'red'))
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


# all the variance in direction of 1st eigenvector?
diag(Lambda.hat)
diag(Lambda.true)[1:J]


image(Sigma.true)
image(Sigma.hat)

par(mfrow=c(2,2))
for (j in 1:J) {
  plot(eigen(Sigma.true)$vectors[,j], Gamma.hat[,j])  
}
par(mfrow=c(1,1))


plot(ts(Y), lwd=2)
lines(-ts(Gamma.hat[,1]*sqrt(diag(Lambda.hat)[1])), col='red', lwd=2)
lines(ts(br.out$fitted.values), lty=2, lwd=2)
legend('topleft', legend=c('True', 'Reg Fit', 'Gamma 1'), col=c('black', 'black', 'red'), lty=c(1,2,1), lwd=2)

plot(Sigma.hat[,1], Sigma.hat[,2])
plot(Y,Sigma.hat[,1])



#plot(ts(Gamma.hat[,1]), ylim=c(-.25,0.25))
#lines(ts(Gamma.hat[,1]+Gamma.hat[,2]), col='green')
#lines(br.out$fitted.values/25, col='red')
#lines(Y/25, col='red', lty=2)

#V.temp <- B.0%*%beta.M2.post%*%t(B.0)
#V.temp.eig <- eigen(V.temp)
#zapsmall(V.temp.eig$values)


#### Freq ####
#m.lm <- lm(Y~B.0-1)
s2 <- 1
beta.cov.freq <- s2*solve(t(B.0)%*%B.0)
beta.hat <- solve(t(B.0)%*%B.0)%*%t(B.0)%*%Y
V.0.freq <- beta.cov.freq #+ beta.hat%*%t(beta.hat)

V.1.freq <- D.0%*%t(Q.0)%*%V.0.freq%*%Q.0%*%D.0
V.1.freq.eig <- eigen(V.1.freq)
P.1.freq <- V.1.freq.eig$vectors
D.1.freq <- diag(V.1.freq.eig$values)


Gamma.hat.freq <- P.0%*%P.1.freq
Lambda.hat.freq <- D.1.freq
Sigma.hat.freq <- Gamma.hat.freq%*%Lambda.hat.freq%*%t(Gamma.hat.freq)

diag(D.1.freq)
image(Sigma.hat.freq)