# Example of using package.
# Data comes from Y ~ GP(0, K(X))




setwd('~/Git/gp_kernlearn/')
source('abc/abc.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')


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
#p <- 1
#X <- matrix(runif(N*p, -5, 5), nrow=N, ncol=p)
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- rbf_kernel(X, 1, 3)
Y <- mvtnorm::rmvnorm(1, rep(0, N), K)

plot(X, Y)
image(K)

X.dm <- make_legendre_design_matrix_1D(3, X)

