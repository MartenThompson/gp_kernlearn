weighted_eval_cov_maker <- function(L, schedule) {
  if (L != length(schedule)) {
    stop('Weight schedule should be same as L (', L, '); is ', length(schedule))
  }
  
  weighted_eval_cov <- function(mat) {
    cov.mat <- cov(mat) 
    e.vals <- eigen(cov.mat, symmetric = TRUE, only.values = TRUE)[[1]][1:L]
    weighted.mean(e.vals, schedule)
  }
  
  return(weighted_eval_cov)
}

max_eigendiscrepancy_cov_maker <- function(obs.mat) {
  e.vals.obs <- eigen(cov(obs.mat), symmetric = TRUE, only.values = TRUE)[[1]]
  
  max_eigendiscrepancy_cov <- function(mat) {
    cov.mat <- cov(mat) 
    e.vals <- eigen(cov.mat, symmetric = TRUE, only.values = TRUE)[[1]]
    return(max(abs(e.vals.obs - e.vals)))
  }
  return(max_eigendiscrepancy_cov)
}

central_difference <- function(values) {
  n <- length(values) 
  cent.diffs <- rep(NA, n) # leave 1 and n NA for clarity
  for (i in 2:(n-1)){
    cent.diffs[i] <- values[i+1] + values[i-1] - 2*values[i]
  }
  return(cent.diffs)
}

max_centraldifflogeigen_cov <- function(mat) {
  cov.mat <- cov(mat)
  e.vals <- eigen(cov.mat, symmetric = TRUE, only.values = TRUE)[[1]]
  return(which.max(central_difference(log(e.vals))))
}


# aka Euclidean norm, sqrt of sum of |a_ij|^2
frob_norm_cov <- function(mat) {
  return(norm(cov(mat), type='F'))
}


# Larges col sum (of abs values)
one_norm_cov <- function(mat) {
  return(norm(cov(mat), type='O'))
}


# Largest row sum (of abs values)
infinity_norm_cov <- function(mat) {
  return(norm(cov(mat), type='I'))
}

# sqrt of max eigval of A^H A
spectral_norm_cov <- function(mat) {
  return(norm(cov(mat), type='2'))
}

eigenvals_cov <- function(mat) {
  e.out <- eigen(cov(mat), symmetric = TRUE, only.values = TRUE)
  return(e.out$values)
}

