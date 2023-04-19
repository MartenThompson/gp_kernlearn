

stan_beta_post <- function(Y, B.0, n.mcmc) {
  br.out <- stan_glm(Y ~ B.0-1, family=gaussian(), iter=n.mcmc)
  beta.post <- matrix(NA, nrow=n.mcmc, ncol=J)
  for (j in 1:J) {
    beta.post[,j] <- unlist(br.out$stanfit@sim$samples[[1]][[j]])
  }
  
  return(list(
    beta_post=beta.post,
    stan_fit=br.out
  ))
}


# Produces n.mcmc estimates of covariance values up to lag (default max possible).
acv_est_single <- function(B.0, beta.post.samples, lag=(N.x-1)) {
  N.x <- dim(B.0)[1]
  
  if (lag >= N.x) {
    stop('Lag cannot be greater than or equal to N.x=', N.x, ', was ', lag)
  }
  
  n.mcmc <- dim(beta.post.samples)[1]
  Y.est <- beta.post.samples %*% t(B.0)
  acv.out <- acf(t(Y.est), plot=FALSE, type='covariance', lag.max = lag)

  # only keep the cov from within a given beta sample (no cross terms)
  acv.ests <- matrix(NA, n.mcmc, lag+1)
  for (i in 1:n.mcmc) {
    acv.ests[i,] <- acv.out$acf[,i,i]
  }
  
  #acv.out GB of likely unneeded info
  return(list(
    acv.ests=acv.ests,
    acv.all=acv.out$acf
  ))
}


est_sigma_acv <- function(Y.list, B.0, n.mcmc, lag=(N.x-1)) {
  N.x <- dim(B.0)[1]
  R <- length(Y.list)
  
  if (lag >= N.x) {
    stop('Lag cannot be greater than or equal to N.x=', N.x, ', was ', lag)
  }
  
  acv.all <- matrix(NA, n.mcmc, lag+1)
  
  for (r in 1:R) {
    cat(r, '/', R, '\n')
    silence.me <- capture.output(
      beta.post.single <- stan_beta_post(Y.list[[r]], B.0, n.mcmc)  
    )
    acv.single <- acv_est_single(B.0, beta.post.single$beta_post, lag)
    
    if (1==r) {
      acv.all[1:n.mcmc,] <- acv.single$acv.ests
    } else {
      acv.all <- rbind(acv.all, acv.single$acv.ests)
    }
  }
  
  return(list(
    acv.all=acv.all
  ))
  
}




make_sigmahat_from_acv <- function(acv.single) {
  n <- length(acv.single)
  Sigma.hat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    if (i < n) {
      idx <- c(i:1,2:(n-i+1))  
    } else {
      idx <- n:1
    }
    
    Sigma.hat[i,] <- acv.single[idx]
  }
  
  return(Sigma.hat)
}