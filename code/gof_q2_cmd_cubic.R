setwd('~/Git/gp_kernlearn/')
library(rstanarm)
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Kern Learn Loop ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

args.in <- commandArgs(trailingOnly = TRUE)
if (2 != length(args.in)){
  stop('Need args: primary degree and R')
}


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
  #image(K.hat)
  
  E.alt <- apply(beta.manysamples.alt, 2, mean)
  V.alt <- cov(beta.manysamples.alt)
  post.ests.alt <- posterior_test_meanvar_brev(E.alt, V.alt, X, Y.history[[3]], X.test, leg_basis_maker_alt)
  K.hat.alt <- post.ests.alt$kern.pieces$K.train
  
  
  return(list(
    K1 = K.hat.pri,
    K2 = K.hat.alt
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Q2 Test ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ni  : number of obs used to est Ki (here, num material samples)
# m   : dimension of response (here, number of x locations)
Q2_stat <- function(K1, K2, N1, N2, m) {
  n1 <- N1 - 1
  n2 <- N2 - 1
  n <- n1 + n2
  V1 <- n1*K1
  V2 <- n2*K2
  V <- V1 + V2
  
  a11 <- (1/(m*n1))*sum(diag(V1))
  a12 <- (1/(m*n2))*sum(diag(V2))
  a21 <- (1/(m*(n1-1)*(n1+2)))*(sum(diag(V1%*%V1)) - (1/n1)*(sum(diag(V1)))^2)
  a22 <- (1/(m*(n2-1)*(n2+2)))*(sum(diag(V2%*%V2)) - (1/n2)*(sum(diag(V2)))^2)
  
  a1 <- (1/(m*n))*sum(diag(V))
  a2 <- (1/(m*(n-1)*(n+2)))*(sum(diag(V%*%V)) - (1/n)*sum(diag(V)^2)) # TODO: Vsquared correct?
  a3 <- (1/(n*(n^2+3*n+4)))*(((1/m)*sum(diag(V%*%V%*%V))) - 3*n*(n+1)*m*a2*a1 - n*m^2*a1^3)
  c0 <- n*(n^3 +6*n^2 + 21*n + 18)
  c1 <- 2*n*(2*n^2 + 6*n +9)
  c2 <- 2*n*(3*n + 2)
  c3 <- n*(2*n^2 + 5*n + 7)
  a4 <- (1/c0)*((1/m)*sum(diag(V%*%V%*%V%*%V)) - m*c1*a1 - m^2*c2*(a1^2)*a2 - m*c3*(a2^2) - n*(m^3)*(a1^4))
  
  gam1 <- a21/(a11^2)
  gam2 <- a22/(a12^2)
  xi1 <- (4/n1^2)*((a2^2/a1^4) + (2*n1/m)*((a2^3/a1^6) - (2*a2*a3/a1^5) + (a4/a1^4)))
  xi2 <- (4/n2^2)*((a2^2/a1^4) + (2*n2/m)*((a2^3/a1^6) - (2*a2*a3/a1^5) + (a4/a1^4)))
  stat <- (gam1 - gam2)/sqrt(xi1 + xi2)
  
  return(list(stat=stat))  
}

#Q2_stat(K, K.temp, n.mat.samp, n.mat.samp, N.x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


N.x <- 10       # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)
K <- cubic_kernel(X, 1, 1, 1, 1, 0.1)
#K <- rbf_kernel(X, 5, 10)
kern.name <- 'cubic1-1-1-1-0.1/' # TODO move all files in quad1-1-1-1-0.1/ into cubic dir

# plot(NA,NA,xlim=c(min(X),max(X)), ylim=c(-8,8))
# for (r in 1:5) {
#   err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), K))
#   mu <- rep(0, N.x)
#   Y <- mu + err
#   points(X, Y, pch=16, col=rgb(runif(1),runif(1),runif(1)))
# }

save_slug <- paste0('package/gp_kernlearn/code/output/gof_q2_', kern.name)
dir.create(file.path(save_slug))

primary.degree <- as.numeric(args.in[1])
alternative.degs <- primary.degree:(primary.degree+1)

n.mat.samp <- as.numeric(args.in[2]) # ie n material samples
n.analysis <- 30 # number of times to get Q2 stat

saveRDS(list(Nx=N.x, X=X,K=K), paste0(save_slug, 'datagen_details_nx', N.x, '_nmat',n.mat.samp,'.RData'))


stat.hist <- matrix(NA, nrow=n.analysis, ncol=length(alternative.degs))
pval.hist <- matrix(NA, nrow=n.analysis, ncol=length(alternative.degs))

cn <- paste0(rep('deg', length(alternative.degs)), alternative.degs)
colnames(stat.hist) <- cn
colnames(pval.hist) <- cn

set.seed(2022)
for (j in 1:length(alternative.degs)) {
  alt.deg <- alternative.degs[j]
  cat('\nAlt Deg', alt.deg, '\n')
  
  for (i in 1:n.analysis) {
    cat(i)
    K.ests <- gen_dat_est_covs(primary.degree, alt.deg, N.x, K, n.mat.samp)
    Q2.out <- Q2_stat(K.ests$K1, K.ests$K2, n.mat.samp, n.mat.samp, N.x)
    
    stat.hist[i,j] <- Q2.out$stat
    pval.hist[i,j] <- 2*pnorm(abs(Q2.out$stat), lower.tail = FALSE)
  }  
}


write.csv(stat.hist, paste0(save_slug, 'stat_hist_prim', primary.degree, '_nx', N.x, '_nmat',n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)
write.csv(pval.hist, paste0(save_slug, 'pval_hist_prim', primary.degree, '_nx', N.x, '_nmat', n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)

out.mat <- matrix(apply(pval.hist, 2, function(col){c(mean(col), signif(mean(col),2), sd(col), signif(sd(col),2))}), nrow=4)
out.mat <- cbind(c('mean', 'mean', 'sd', 'sd'), out.mat)
colnames(out.mat) <- c('', cn)
write.csv(out.mat, paste0(save_slug, 'pval_sum_prim', primary.degree, '_nx', N.x, '_nmat', n.mat.samp, '_nsamp', n.analysis, '.csv'), row.names = FALSE)