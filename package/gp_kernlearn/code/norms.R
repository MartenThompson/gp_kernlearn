

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
