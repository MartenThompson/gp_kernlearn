require(MonoPoly)

dist_euc <- function(vect) {
  D <- as.matrix(dist(vect, method = 'euclidean', diag = TRUE))
  return(D[upper.tri(D, diag = TRUE)])
}

modfit_lm <- function(Y, X) {
  return(lm(Y~X))
}

modfit_lm3 <- function(Y, X) {
  return(lm(Y ~ X + I(X^2) + I(X^3)))
}

modfit_monopoly_maker <- function(degree) {
  modfit_monopoly <- function(Y, X) {
    return(monpol(Y ~ X, degree=degree))
  }
  
  return(modfit_monopoly)
}


#modfit_spline <- function(Y, X) {
#  return(smooth.spline(x=X, y=Y, nknots = 1))
#}

# reg_fitter takes args (Y, X)
make_stationary <- function(K, X, dist_metric, reg_fitter) {
  K.stationary <- matrix(NA, nrow=dim(K)[1], ncol=dim(K)[2])
  K.uptri <- K[upper.tri(K, diag=TRUE)]
  D <- dist_metric(X) 
  
  m <- reg_fitter(K.uptri, D) 
  K.stat.pred <- predict(m, data.frame(X=D))
  
  K.stationary[upper.tri(K.stationary, diag=TRUE)] <- K.stat.pred
  K.stationary[lower.tri(K.stationary)] = t(K.stationary)[lower.tri(K.stationary)]
  
  if(!isSymmetric(K.stationary)) {
    stop('K.stationary not symmetric')
    return(K.stationary)
  }
  
  return(list(
    K=K, dist_metric=dist_metric, reg_fitter=reg_fitter,
    D=D, model=m, modtrain=K.uptri, modpred=K.stat.pred, 
    K.stationary=K.stationary
  ))
}

