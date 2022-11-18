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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### p = 1. Y=poly(x) + GP(0, rbf) rough ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
save_slug <- 'package/gp_kernlearn/vignettes/output/br_polymean_rbf_p1/'
dir.create(file.path(save_slug))

N <- 50
set.seed(2022)
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- rbf_kernel(X, 1, 5, 0.1)
err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
mu <- 0*(-2*X + 0.5*X^2 + 0.1*X^3) 
Y <- mu + err
#saveRDS(data.frame(X=X,Y=Y,mu=mu,err=err), paste0(save_slug, 'data.RData'))
#saveRDS(K, paste0(save_slug, 'dataK.RData'))

image(K)
plot(X, Y)
points(X, mu, col='green')
points(X, err, col='red')

K[1:5,1:5]

leg.deg <- 4
leg_basis_maker <- make_legendre1D_basis_maker(degree=leg.deg)
b.X <- leg_basis_maker(X)

br.output <- stan_glm(Y ~ b.X[,2:(leg.deg+1)], family=gaussian())
print(br.output)

# correctly identified when sigma=0?!
#saveRDS(br.output, paste0(save_slug, 'br_output_10min_deg3.RData'))

#beta.samples <- matrix(unlist(br.output$stanfit@sim$samples[[1]]), ncol=length(br.output$stanfit@sim$samples[[1]]))

beta.samples <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[1]]), ncol=1)
for (i in 2:(leg.deg+1)) {
  new <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[i]]), ncol=1)
  beta.samples <- cbind(beta.samples, new)
}

cov(beta.samples)
br.output$covmat

apply(beta.samples, 2, mean)
br.output$coefficients

par(mfrow=c(2,3))
for (i in 1:(leg.deg+1)) {
  plot(density(br.output$stanfit@sim$samples[[1]][[i]]), main=paste0('Beta ', i))
  #lines(c(E.leg[i], E.leg[i]), c(0,1), col='red')
}
par(mfrow=c(1,1))


X.test <- matrix(sort(runif(100, -5, 5)), ncol=1)
post.ests <- posterior_test_meanvar_br(br.output, X, Y, X.test, leg_basis_maker)

plot(X, Y)
points(X, mu, col='green')
points(X.test, post.ests$post.test.mean, col='blue', pch=16)
points(rep(X.test,2), c(post.ests$post.test.mean + sqrt(diag(post.ests$kern.pieces$K.star)) + 0.9,
                        post.ests$post.test.mean - sqrt(diag(post.ests$kern.pieces$K.star)) - 0.9),
       pch=16, col=rgb(0,0,1,0.5))
image(post.ests$kern.pieces$K.train)
image(post.ests$kern.pieces$K.star)

K[1:6,1:6]
post.ests$kern.pieces$K.train[1:6,1:6]
image(K)
image(post.ests$kern.pieces$K.train + diag(0.5, nrow=dim(K)[1], ncol=dim(K)[1]))
#TODO calc matrix norms of K and K learned (w/ sigma^2) over differing leg.degree


image(K - post.ests$kern.pieces$K.train)

frobenius_norm <- function(mat) {
  sum(mat^2)
}

frobenius_norm(K - post.ests$kern.pieces$K.train)
# deg1 616
# deg3 620
# deg7 625
# deg10 625




