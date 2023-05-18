library(rstanarm)
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

args.in <- commandArgs(trailingOnly = TRUE)
if (2 != length(args.in)){
  stop('Need args: primary degree and R')
}



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

primary.degree <- as.numeric(args.in[1])
alternative.degs <- primary.degree:(primary.degree+1)

n.mat.samp <- as.numeric(args.in[2]) # ie n material samples
n.analysis <- 30 # number of times to get box stat

saveRDS(list(Nx=N.x, X=X,K=K), paste0(save_slug, 'datagen_details_nx', N.x, '_nmat',n.mat.samp,'.RData'))

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