# Does Box's M Test mesh with our method?
# Can we reliably test Sigma1 = Sigma2 as a proxy for 
# basis degree #1 being sufficient to express data.


# IDK how this applies to higher basis dimension, but...
# Generate data according to GP(0, lin_kern).
# We want to investigate behavior of Box's M stat as we compare 
# deg1 = 0 to deg2 = 0,1,2,3,4,5.


setwd('~/Git/gp_kernlearn/')
library(rstanarm)
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Kern Learn Loop ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# assumes 0 centered X
lin_kernel <- function(X, sigma.0, sigma.1, fuzz=0.1) {
  K <- sigma.0 + sigma.1*X%*%t(X)
  diag(K) <- diag(K) + fuzz
  return(K)
}


gen_dat_est_covs <- function(leg.deg, reps){
  leg_basis_maker <- make_legendre1D_basis_maker(degree=leg.deg)
  b.X <- leg_basis_maker(X)
  
  leg_basis_maker_deg1 <- make_legendre1D_basis_maker(degree = 1)
  b1.X <- leg_basis_maker_deg1(X)
  
  beta.samples.deg1 <- NA
  beta.samples.deg.alt <- NA
  Y.history <- list()
  err.history <- list()

  for (r in 1:reps) {
    #cat('Replication', r, '\n')
    err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
    mu <- rep(0, N)
    Y <- mu + err
    Y.history[[r]] <- Y
    err.history[[r]] <- err
    
    silent <- capture.output(br.output <- stan_glm(Y ~ b.X[,2:(leg.deg+1)], family=gaussian()))
    beta.samples <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[1]]), ncol=1)
    for (i in 2:(leg.deg+1)) {
      new <- matrix(unlist(br.output$stanfit@sim$samples[[1]][[i]]), ncol=1)
      beta.samples <- cbind(beta.samples, new)
    }
    
    silent <- capture.output(br1.output <- stan_glm(Y ~ b1.X[,2], family = gaussian()))
    beta1.samples <-  matrix(unlist(br1.output$stanfit@sim$samples[[1]][[1]]), ncol=1)
    new <- matrix(unlist(br1.output$stanfit@sim$samples[[1]][[2]]), ncol=1)
    beta1.samples <- cbind(beta1.samples, new)
    
    if (1==r) {
      beta.samples.deg.alt <- beta.samples
      beta.samples.deg1 <- beta1.samples
    } else {
      beta.samples.deg.alt <- rbind(beta.samples.deg.alt, beta.samples)
      beta.samples.deg1 <- rbind(beta.samples.deg1, beta1.samples)
    }
  }
  
  #saveRDS(Y.history, paste0(save_slug, 'leg', leg.deg, '/Yhist.RData'))
  #saveRDS(err.history, paste0(save_slug, 'leg', leg.deg, '/errhist.RData'))
  #saveRDS(beta.samples.many, paste0(save_slug, 'leg', leg.deg, '/betasamples.RData'))
  
  E <- apply(beta.samples.deg.alt, 2, mean)
  V <- cov(beta.samples.deg.alt)
  #saveRDS(E, paste0(save_slug, 'leg', leg.deg, '/E.RData'))
  #saveRDS(V, paste0(save_slug, 'leg', leg.deg, '/V.RData'))
  
  X.test <- matrix(seq(-5,5,length.out=20), ncol=1)
  post.ests <- posterior_test_meanvar_brev(E, V, X, Y.history[[3]], X.test, leg_basis_maker)
  K.hat <- post.ests$kern.pieces$K.train
  
  E.deg1 <- apply(beta.samples.deg1, 2, mean)
  V.deg1 <- cov(beta.samples.deg1)
  post.ests.deg1 <- posterior_test_meanvar_brev(E.deg1, V.deg1, X, Y.history[[3]], X.test, leg_basis_maker_deg1)
  K.hat.deg1 <- post.ests.deg1$kern.pieces$K.train
  
  
  return(list(
    K1 = K.hat,
    K2 = K.hat.deg1
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Box's M Test ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reps: n_j
# n.x : k
boxM_stat <- function(K1, K2, reps, n.x) {
  n <- 2*reps
  S.pooled = (1/(n-2))*(((reps-1)*K1) + ((reps-2)*K2))  
  M = (n-2)*log(det(S.pooled)) - (reps-1)*log(det(K1)) - (reps-1)*log(det(K2))
  c = ((2*n.x^2 + 3*n.x -1)/(6*(n.x+1)))*((2/(reps-1)) - (1/(n-2)))
  stat = M*(1-c)
  df = 0.5*n.x*(n.x+1)
  
  return(list(stat=stat,
              df = df,
              n=n,
              S.pooled = S.pooled,
              M = M,
              c = c))  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save_slug <- 'package/gp_kernlearn/code/output/gof_boxm/'
dir.create(file.path(save_slug))


N <- 10
X <- matrix(seq(-5,5,length.out=N), nrow=N, ncol=1)
K <- lin_kernel(X, 1, 1/25, 0.1)

n.mat.samp <- 20 # e.g. n material samples
n.analysis <- 1 # number of times to get box stat

leg.degs <- 1:2

stat.hist <- matrix(NA, nrow=n.analysis, ncol=length(leg.degs))
df.hist <- matrix(NA, nrow=n.analysis, ncol=length(leg.degs))
pval.hist <- matrix(NA, nrow=n.analysis, ncol=length(leg.degs))

colnames(stat.hist) <- paste0(rep('deg', length(leg.degs)), leg.degs)
colnames(pval.hist) <- paste0(rep('deg', length(leg.degs)), leg.degs)
colnames(df.hist) <- paste0(rep('deg', length(leg.degs)), leg.degs)

set.seed(2022)
for (j in 1:length(leg.degs)) {
  leg.deg <- leg.degs[j]
  cat('\nLeg Deg', leg.deg, '\n')
  
  for (i in 1:n.analysis) {
    cat(i)
    K.ests <- gen_dat_est_covs(leg.deg, n.mat.samp)
    boxM.out <- boxM_stat(K.ests$K1, K.ests$K2, n.mat.samp, N)
    
    stat.hist[i,j] <- boxM.out$stat
    df.hist[i,j] <- boxM.out$df
    pval.hist[i,j] <- pchisq(boxM.out$stat, boxM.out$df, lower.tail = FALSE)
  }  
}


write.csv(stat.hist, paste0(save_slug, 'stat_hist_nmat',n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)
write.csv(df.hist, paste0(save_slug, 'df_hist_nmat',n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)
write.csv(pval.hist, paste0(save_slug, 'pval_hist_nmat',n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)

