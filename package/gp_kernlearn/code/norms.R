

frob_norm_cov <- function(mat) {
  return(norm(cov(mat), type='F'))
}

one_norm_cov <- function(mat) {
  return(norm(cov(mat), type='O'))
}

infinity_norm_cov <- function(mat) {
  return(norm(cov(mat), type='I'))
}

spectral_norm_cov <- function(mat) {
  return(norm(cov(mat), type='2'))
}