# TODO: normalizing needs to be done in a global fashion. Currently just a few datasets informed scaling factors.
read_normalize_unit_square <- function(filename) {
  raw <- read.csv(paste0('ms_eg/jan2022_bcc_files_trimmed/',filename), sep=' ')
  Y <- raw$G.J.mol.-raw$G_MIX.J.mol.
  Y <- (Y + 179569.5)/124340.5
  T <- (raw$T.K.- 1920)/1080
  C <- (raw$X.A. - 0.5)*2
  
  return(list(Y=Y,
              X=matrix(c(T,C), ncol=2)))
}