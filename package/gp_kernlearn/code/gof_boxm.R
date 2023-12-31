# Does Box's M Test mesh with our method?
# Can we reliably test Sigma1 = Sigma2 as a proxy for 
# basis degree #1 being sufficient to express data.


# IDK how this applies to higher basis dimension, but...
# Generate data according to GP(0, lin_kern).
# We want to investigate behavior of Box's M stat as we compare 
# deg1 = 1 to deg2 = 1,2,3,4,5.
# If we fail to reject equality, there's no need to consider higher deg.

# Results
# For linear kernel data
  # failed to reject deg1(1) == deg2(1); barel rej d1(1) == d2(2); rej d1(1) == d2(3). So 1 is slightly diff than 2, v diff than 3.   
  # not rej d1(2) == d2(2); barely rej d1(2)==d2(3); rej d1(2)==d2(4).
# For quad kernel data
  # failed to reject deg1(1) == deg2(1); rejected deg1(1) == deg2(2, 3, 4, 5). So 1 is different than 2,3,4,5. 
  # failed to reject deg1(2) == deg2(2); barely rej d1(2)==d2(3); rej d1(2)==d2(4). Maybe something here?

# Not satisfactory. What we want to see is:
# For linear data: 1 not that different than 2,3 so don't bother.
# For quad data: 1 different than 2; 2 not that diff than 3, 4 so go with 2.

rm(list=ls())

setwd('~/Git/gp_kernlearn/')
library(rstanarm)
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Kern Learn Loop ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gen_dat_est_covs <- function(pri.deg, alt.deg, N, K, reps){
  leg_basis_maker_pri <- make_legendre1D_basis_maker(degree = pri.deg)
  b.X.pri <- leg_basis_maker_pri(X)
  
  leg_basis_maker_alt <- make_legendre1D_basis_maker(degree=alt.deg)
  b.X.alt <- leg_basis_maker_alt(X)
  
  beta.manysamples.pri <- NA
  beta.manysamples.alt <- NA
  Y.history <- list()
  err.history <- list()
  
  for (r in 1:reps) {
    err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
    mu <- rep(0, N)
    Y <- mu + err
    Y.history[[r]] <- Y
    err.history[[r]] <- err
    
    silent <- capture.output(br.out.pri <- stan_glm(Y ~ b.X.pri[,2:(pri.deg+1)], family=gaussian()))
    beta.samples.pri <- matrix(unlist(br.out.pri$stanfit@sim$samples[[1]][[1]]), ncol=1)
    for (i in 2:(pri.deg+1)) {
      new <- matrix(unlist(br.out.pri$stanfit@sim$samples[[1]][[i]]), ncol=1)
      beta.samples.pri <- cbind(beta.samples.pri, new)
    }
    
    silent <- capture.output(br.out.alt <- stan_glm(Y ~ b.X.alt[,2:(alt.deg+1)], family=gaussian()))
    beta.samples.alt <- matrix(unlist(br.out.alt$stanfit@sim$samples[[1]][[1]]), ncol=1)
    for (i in 2:(alt.deg+1)) {
      new <- matrix(unlist(br.out.alt$stanfit@sim$samples[[1]][[i]]), ncol=1)
      beta.samples.alt <- cbind(beta.samples.alt, new)
    }
    
    if (1==r) {
      beta.manysamples.pri <- beta.samples.pri
      beta.manysamples.alt <- beta.samples.alt
    } else {
      beta.manysamples.pri <- rbind(beta.manysamples.pri, beta.samples.pri)
      beta.manysamples.alt <- rbind(beta.manysamples.alt, beta.samples.alt)
    }
  }
  
  E.pri <- apply(beta.manysamples.pri, 2, mean)
  V.pri <- cov(beta.manysamples.pri)
  
  X.test <- matrix(seq(-5,5,length.out=20), ncol=1)
  post.ests.pri <- posterior_test_meanvar_brev(E.pri, V.pri, X, Y.history[[3]], X.test, leg_basis_maker_pri)
  K.hat.pri <- post.ests.pri$kern.pieces$K.train
  
  E.alt <- apply(beta.manysamples.alt, 2, mean)
  V.alt <- cov(beta.manysamples.alt)
  post.ests.alt <- posterior_test_meanvar_brev(E.alt, V.alt, X, Y.history[[3]], X.test, leg_basis_maker_alt)
  K.hat.alt <- post.ests.alt$kern.pieces$K.train
  # image(K.hat.alt)
  
  return(list(
    K1 = K.hat.pri,
    K2 = K.hat.alt
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Box's M Test ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://www.real-statistics.com/multivariate-statistics/boxs-test/boxs-test-basic-concepts/
# https://rdrr.io/cran/rstatix/src/R/box_m.R
# reps: n_j
# n.x : k
boxM_stat <- function(K1, K2, reps, n.x) {
  n <- 2*reps
  S.pooled = (1/(n-2))*(((reps-1)*K1) + ((reps-1)*K2))  
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
N.x <- 10         # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)
K <- lin_kernel(X, 1, 1, 0.1)
#K <- quad_kernel(X, 1, 1/25, 1/25, 0.1)
#K <- cubic_kernel(X, 1, 1/25, 1/25, 1/25, 0.1)
kern.name <- 'lin1-1-0.1/'

save_slug <- paste0('package/gp_kernlearn/code/output/gof_boxm_', kern.name)
dir.create(file.path(save_slug))


plot(NA,NA,xlim=c(min(X),max(X)), ylim=c(-8,8))
for (r in 1:5) {
  err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), K))
  mu <- rep(0, N.x)
  Y <- mu + err
  points(X, Y, pch=16, col=rgb(runif(1),runif(1),runif(1)))
}


n.mat.samp <- 20 # e.g. n material samples
n.analysis <- 2 # number of times to get box stat

saveRDS(list(Nx=N.x, X=X,K=K), paste0(save_slug, 'datagen_details_nx', N.x, '_nmat',n.mat.samp,'.RData'))

primary.degree <- 1
alternative.degs <- 1:2

stat.hist <- matrix(NA, nrow=n.analysis, ncol=length(alternative.degs))
df.hist <- matrix(NA, nrow=n.analysis, ncol=length(alternative.degs))
pval.hist <- matrix(NA, nrow=n.analysis, ncol=length(alternative.degs))

cn <- paste0(rep('deg', length(alternative.degs)), alternative.degs)
colnames(stat.hist) <- cn
colnames(pval.hist) <- cn
colnames(df.hist) <- cn

set.seed(2022)
for (j in 1:length(alternative.degs)) {
  alt.deg <- alternative.degs[j]
  cat('\nAlt Deg', alt.deg, '\n')
  
  for (i in 1:n.analysis) {
    cat(i)
    K.ests <- gen_dat_est_covs(primary.degree, alt.deg, N.x, K, n.mat.samp)
    boxM.out <- boxM_stat(K.ests$K1, K.ests$K2, n.mat.samp, N.x)
    
    stat.hist[i,j] <- boxM.out$stat
    df.hist[i,j] <- boxM.out$df
    pval.hist[i,j] <- pchisq(boxM.out$stat, boxM.out$df, lower.tail = FALSE)
  }  
}


write.csv(stat.hist, paste0(save_slug, 'stat_hist_prim', primary.degree, '_nx', N.x, '_nmat',n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)
write.csv(df.hist, paste0(save_slug, 'df_hist_prim', primary.degree, '_nx', N.x, '_nmat', n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)
write.csv(pval.hist, paste0(save_slug, 'pval_hist_prim', primary.degree, '_nx', N.x, '_nmat', n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)

out.mat <- matrix(apply(pval.hist, 2, function(col){c(mean(col), signif(mean(col),2), sd(col), signif(sd(col),2))}), nrow=4)
out.mat <- cbind(c('mean', 'mean', 'sd', 'sd'), out.mat)
colnames(out.mat) <- c('', cn)
write.csv(out.mat, paste0(save_slug, 'pval_sum_prim', primary.degree, '_nx', N.x, '_nmat', n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)