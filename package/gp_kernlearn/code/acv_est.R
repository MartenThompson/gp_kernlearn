

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
acv_est_single <- function(B.0, beta.post.samples, lag=(N.x-1), type='covariance') {
  N.x <- dim(B.0)[1]
  
  if (lag >= N.x) {
    stop('Lag cannot be greater than or equal to N.x=', N.x, ', was ', lag)
  }
  
  n.mcmc <- dim(beta.post.samples)[1]
  Y.est <- beta.post.samples %*% t(B.0)
  
  acv.ests <- t(apply(Y.est, 1, function(row) {acf(row, plot=FALSE, type=type, lag.max=lag)$acf}))
  
  #full output of acf is GB of likely unneeded info
  return(list(
    acv.ests=acv.ests,
    Y.est = Y.est
  ))
}


est_acv_within <- function(Y.list, B.0, n.mcmc, lag=(N.x-1), type='covariance') {
  N.x <- dim(B.0)[1]
  R <- length(Y.list)
  
  if (lag >= N.x) {
    stop('Lag cannot be greater than or equal to N.x=', N.x, ', was ', lag)
  }
  
  acv.all <- matrix(NA, n.mcmc, lag+1) # ends up being n.mcmc*R rows
  yhat.all <- matrix(NA, n.mcmc, N.x)  # ends up being n.mcmc*R rows
  
  # Get within sample (r) variance for all n.mcmc-many predicted y.hats.
  for (r in 1:R) {
    cat(r, '/', R, '\n')
    silence.me <- capture.output(
      beta.post.single <- stan_beta_post(Y.list[[r]], B.0, n.mcmc)  
    )
    
    acv.single <- acv_est_single(B.0, beta.post.single$beta_post, lag, type)
    resids <- t(apply(acv.single$Y.est, 1, function(row){row-Y.list[[r]]}))
    
    if (1==r) {
      acv.all[1:n.mcmc,] <- acv.single$acv.ests
      yhat.all[1:n.mcmc,] <- acv.single$Y.est
    } else {
      acv.all <- rbind(acv.all, acv.single$acv.ests)
      yhat.all <- rbind(yhat.all, acv.single$Y.est)
    }
  }
  
  # Get between sample variance, having n.mcmc-many groups each with R values.
  yhat.all.means <- apply(yhat.all, 1, mean)
  acf.bt.mcmc <- matrix(NA, nrow=n.mcmc, ncol=(lag+1))
  for (i in 1:n.mcmc) {
    idx.group <- seq(i, R*n.mcmc, by=n.mcmc)
    y.means.group <- yhat.all.means[idx.group]
    acf.bt.mcmc[i,] <- acf(y.means.group, lag.max = lag, type='covariance', plot=FALSE)$acf[,1,1]
  }
  
  # Combine within and between sample variance. 
  # We have R*n.mcmc 'within' acf vectors (ie acv.all), and we have n.mcmc 'between'
  # acf vectors (ie acf.bt.mcmc).
  acv.all.w.bt <- acv.all
  for (i in 0:(nrow(acv.all)-1)) {
    acv.all.w.bt[i+1,] <- acv.all.w.bt[i+1,]+acf.bt.mcmc[(i%%n.mcmc)+1]
  }
  
  return(list(
    acv.all=acv.all,
    acv.all.w.bt=acv.all.w.bt,
    yhat.all=yhat.all,
    resids=resids
  ))
  
}



est_acv_bt <- function(Y.list, B.0, n.mcmc, lag=(N.x-1), type='covariance') {
  N.x <- dim(B.0)[1]
  R <- length(Y.list)
  
  if (lag >= N.x) {
    stop('Lag cannot be greater than or equal to N.x=', N.x, ', was ', lag)
  }
  
  Y.hat.all <- matrix(NA, nrow=n.mcmc, ncol=N.x) # eventually R*N.x cols
  
  for (r in 1:R) {
    cat(r, '/', R, '\n')
    silence.me <- capture.output(
      beta.post.single <- stan_beta_post(Y.list[[r]], B.0, n.mcmc)  
    )
    
    if (1==r) {
      Y.hat.all[1:n.mcmc,1:N.x] <- beta.post.single$beta_post %*% t(B.0)
    } else {
      Y.hat.all <- cbind(Y.hat.all, beta.post.single$beta_post %*% t(B.0))
    }
  }
  
  acv.ests <- t(apply(Y.hat.all, 1, function(row) {acf(row, plot=FALSE, type=type, lag.max=lag)$acf}))
  
  return(list(
    acv.all=acv.ests
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



#(var.t <- var(unlist(Y.train)))
#var.w <- sum(sapply(1:R, function(r){(1/R)*var(Y.train[[r]])}))
#var.b <- sum(sapply(1:R, function(r){(1/R)*(mean(Y.train[[r]]) - mean(unlist(Y.train)))^2}))
#var.w+var.b


S_B <- function(Y.all, r, s) {
  T <- length(unlist(Y.all)) 
  Y.T.mean <- mean(unlist(Y.all))
  Y.r.mean <- mean(Y.all[[r]])
  Y.s.mean <- mean(Y.all[[s]])
  Y.t <- unlist(c(Y.all[[r]], Y.all[[s]]))
  ss <- sum(Y.t-Y.T.mean)
  #(1/T)*
  return((Y.r.mean-Y.s.mean)*(ss))
}

#S_B(Y.train, 1,4)
#lag <- 5
#sum(sapply((1+lag):R, function(r){S_B(Y.train,r, r-lag)}))
#sum(sapply((1+lag):R, function(r){S_B(Y.train,r-lag, r)}))



#Y.ea.mean <- sapply(1:R, function(r){mean(Y.list[[r]])})
#acf.bt <- acf(Y.ea.mean, lag.max=lag, type=type, plot=FALSE)$acf[,1,1]