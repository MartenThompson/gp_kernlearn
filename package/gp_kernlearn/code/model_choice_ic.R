mvn_loglikelihood <- function(Y, mu, Sigma) {
  mvtnorm::dmvnorm(Y, mean=mu, sigma=Sigma, log=TRUE)
}

aic <- function(llik, n.params) {
  return(-2*llik + 2*n.params)
}

bic <- function(llik, n.params, n.obs) {
  return(-2*llik + log(n.obs)*n.params)
}


