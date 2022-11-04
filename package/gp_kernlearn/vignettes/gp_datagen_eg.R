# Example of using package.
# Data comes from Y ~ GP(0, K(X))




setwd('~/Git/gp_kernlearn/')
source('abc/abc.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

rbf_kernel <- function(X, alpha.0, alpha.1, phi=0) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  if (p >= N) {
    stop("Have not thought through p > N.")
  }
  
  K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
  diag(K) <- diag(K) + phi
  return(K)
}


N <- 50

set.seed(2022)
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- rbf_kernel(X, 1, 5)
Y <- t(mvtnorm::rmvnorm(1, rep(0, N), K))

plot(X, Y)
image(K)

leg_basis_maker <- make_legendre1D_basis_maker(degree=4)
b.X <- leg_basis_maker(X)
#beta.true <- c(-1,3,1,3)
data_gen <- data_gen_maker(b.X)
# TODO: prior hard coded 5
garbage <- data_gen(NA, param_prior()) # just need to run this method once o/w lazy init bites us.
#Y <- data_gen(NA, list(beta=beta.true)) 
#plot(X, Y)
 
M <- 100
runtime <- 5
abc.output <- abc_knn_fixedrt(M, runtime, param_prior, data_gen, Y, distance, n_cores=6, k=1, packages_load=c())

par(mfrow=c(2,3))
for (i in 1:5) {
 plot(density(abc.output$results[,i]), main=paste0('Beta ', i))
 #lines(c(E.leg[i], E.leg[i]), c(0,1), col='red')
}
par(mfrow=c(1,1))

# TODO move to kern
V <- cov(as.matrix(abc.output$results[,1:5])) # TODO: 5 hard coded
E <- apply(abc.output$results[,1:5], 2, mean) # TODO
X.test <- matrix(runif(25, -5, 5), ncol=1)
kern.pieces <- make_kern(X, X.test, V, leg_basis_maker)
post.ests <- posterior_test_meanvar(X.test, X, kern.pieces$K.traininv, 
                                    kern.pieces$K.startr, kern.pieces$K.star,
                                    E, Y, leg_basis_maker)

plot(X, Y)
points(X.test, post.ests$post.test.mean, col='green', pch=16)
image(kern.pieces$K.train)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### p > 1 Example ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#p <- 1
#X <- matrix(runif(N*p, -5, 5), nrow=N, ncol=p)



