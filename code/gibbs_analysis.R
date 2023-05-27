rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/R/kernlearn.R')
source('package/gp_kernlearn/R/basis_orthog_poly.R')
source('package/gp_kernlearn/R/acv_est.R')
source('code/gibbs_vis.R')

require(plotly)

skinny <- function(kernlearn_out) {
  temp <- kernlearn_out
  temp$stan.output <- kernlearn_out$stan.output$stanfit
  return(temp)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Import ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
root.path <- './code/output/gibbs/'
Y.all <- readRDS(paste0('./code/output/gibbs/', 'Y_scaled.Rdata'))
X.all <- readRDS(paste0('./code/output/gibbs/', 'X_scaled.Rdata'))

leg.deg <- 2
basis_maker <- make_legendre2D_basis_maker(leg.deg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Train ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.each <- as.numeric(lapply(Y.all, length))

Y.all.ord <- list()
X.all.ord <- list()
rep.itr <- 1
for (i in 1:length(n.each)) {
  idx.small <- which(n.each == sort(n.each)[i])
  if (1 < length(idx.small)) {
    n.equal.length <- length(idx.small)
    idx.small <- idx.small[rep.itr]
    
    if (rep.itr == n.equal.length) {
      rep.itr <- 1
    } else {
      rep.itr <- rep.itr + 1  
    }
    
  } else {
    rep.itr <- 1
  }
  cat(idx.small, ' ')  
  Y.all.ord[[i]] <- Y.all[[idx.small]]
  X.all.ord[[i]] <- X.all[[idx.small]]
}
Y.all.ord <- Y.all.ord[-1] # too small, not interesting
X.all.ord <- X.all.ord[-1]
plot(ts(as.numeric(lapply(Y.all.ord, length))))

n.train <- 75 #floor(0.5*n.total)
b.X.train <- list()
X.train.long <- X.all.ord[[1]]
Y.train.long <- Y.all.ord[[1]]
for (i in 1:n.train) {
  cat(i,' ')
  b.X.train[[i]] <- basis_maker(X.all.ord[[i]])
  
  if (1 != i) {
    X.train.long <- rbind(X.train.long, X.all.ord[[i]])
    Y.train.long <- c(Y.train.long, Y.all.ord[[i]])
  }
}

stan.chains <- 2
stan.iters <- 200
kernlearn.out <- fit_kernlearn_Xunique(b.X.train, basis_maker, Y.all.ord[1:n.train], X.train.long, Y.train.long, 
                                       stan.chains=stan.chains, stan.iter=stan.iters, EV.only = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Save ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.path <- paste0(root.path, 'analysis_deg', leg.deg, '_ntr', n.train, '_stan', stan.chains, '-', stan.iters, '/')
dir.create(save.path)

file.copy('./package/gp_kernlearn/code/gibbs_analysis.R', paste0(save.path, 'gibbs_analysis.R'))

plot3D_many(X.all.ord[1:n.train], Y.all.ord[1:n.train])

saveRDS(list(leg.deg=leg.deg,
             n.train=n.train,
             b.X.train=b.X.train,
             X.train.long=X.train.long,
             Y.train.long=Y.train.long,
             skinny(kernlearn.out))
        , paste0(save.path, 'analysis_ntr', n.train, '_stan', stan.chains, '-', stan.iters, '.Rdata'))

E.beta <- kernlearn.out$E.est
V.beta <- kernlearn.out$V.est

png(paste0(save.path, 'Vbeta.png'), height=5, width=5, units='in', res=100)
image(V.beta[nrow(V.beta):1,], xaxt='n', yaxt='n')
dev.off()

test.itr <- 100 #n.train + 1
X.test <- X.all[[test.itr]]
Y.test <- Y.all[[test.itr]]


summary(X.test)
x.temp.q50 <- quantile(X.test[,1], 0.5)
X.test.templow <- X.test[X.test[,1] <= x.temp.q50, ]
X.test.temphigh <- X.test[X.test[,1] > x.temp.q50, ]
Y.test.templow <- Y.test[X.test[,1] <= x.temp.q50]
Y.test.temphigh <- Y.test[X.test[,1] > x.temp.q50]


#post.out.obs.templow.skinnycond <- posterior_test_meanvar_brev(E.beta, V.beta,
#                                                               X.test.templow, Y.test.templow,
#                                                               X.test.temphigh, basis_maker)

#plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow.skinnycond$post.test.mean, X.test.templow, Y.test.templow)

X.test.grid <- as.matrix(expand.grid(seq(min(X.test[,1]),max(X.test[,1]), length.out=30),
                                     seq(min(X.test[,2]),max(X.test[,2]), length.out=30), KEEP.OUT.ATTRS = FALSE))
post.out.obs.templow.grid <- posterior_test_meanvar_brev(E.beta, V.beta,
                                                         X.test.templow, Y.test.templow,
                                                         X.test.grid, basis_maker)
plot3D_many(list(X.test.templow,X.test.temphigh, X.test.grid), 
            list(Y.test.templow, Y.test.temphigh, post.out.obs.templow.grid$post.test.mean), 
            c('Conditioned', 'Held Out', 'Post. Mean'), 
            c(I('gray'), I('black'), I('green')))


