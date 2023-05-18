rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
save.slug <- 'package/gp_kernlearn/code/output/data_gen_vis/'


dgobj <- function(name, legend, K) {
  return(list(name=name, 
              legend=legend,
              K=K))
}

plot_nice <- function(n.viz, datagen.list, ylim=c(-5,5)) {
  plot.colors <- viridis::magma(length(datagen.list)+2,0.75)[2:(2+length(datagen.list))]
  
  plot(NA, NA, xlim=c(-1,1), ylim=ylim, ylab=NA, xlab='X')
  for (i in 1:length(datagen.list)) {
    for (r in 1:n.viz) {
      Y <- mvtnorm::rmvnorm(1, sigma=datagen.list[[i]]$K)
      lines(X, Y, lwd=3, col=plot.colors[i])  
    }
  }
  legend('topright', legend=lapply(datagen.list, function(l){l$legend}), col = plot.colors, lty=1, lwd=3)  
}

N.x <- 100         # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)

datagen.list <- list()
datagen.list[[1]] <- dgobj('lin1-1-fz0', '(1,1)', lin_kernel(X, 1, 1, 0))
datagen.list[[2]] <- dgobj('lin5-1-fz0', '(5,1)', lin_kernel(X, 5, 1, 0))
datagen.list[[3]] <- dgobj('lin1-5-fz0', '(1,5)', lin_kernel(X, 1, 5, 0))
datagen.list[[4]] <- dgobj('lin5-5-fz0', '(5,5)', lin_kernel(X, 5, 5, 0))

plot_nice(3, datagen.list)

datagen.list <- list()
datagen.list[[1]] <- dgobj('quad1-1-1-fz0', '(1,1,1)', quad_kernel(X,1, 1, 1, 0))
datagen.list[[2]] <- dgobj('quad5-1-1-fz0', '(5,1,1)', quad_kernel(X,5, 1, 1, 0))
datagen.list[[3]] <- dgobj('quad1-5-1-fz0', '(1.5,1)', quad_kernel(X,1,5,  1, 0))
datagen.list[[4]] <- dgobj('quad1-1-5-fz0', '(1,1,5)', quad_kernel(X,1, 1, 5, 0))
plot_nice(3, datagen.list)

datagen.list <- list()
datagen.list[[1]] <- dgobj('lin1-5-fz0', 'linear', lin_kernel(X, 1, 1, 0.1))
datagen.list[[2]] <- dgobj('quad1-1-5-fz0', 'quadratic', quad_kernel(X,1, 1, 1, 0.1))
datagen.list[[3]] <- dgobj('cubic1-1-1-5-fz0', 'cubic', cubic_kernel(X, 1, 1, 1, 1, 0.1))
png(paste0(save.slug, 'linquadcubic_fz01_dens.png'), height=6, width=6, units='in', res = 100)
plot_nice(2, datagen.list, c(-3,3))
dev.off()



datagen.list <- list()
datagen.list[[1]] <- dgobj('rbf1-0.05', '(1,0.05)', rbf_kernel(X, 1, 0.05, 0.1))
datagen.list[[2]] <- dgobj('rbf5-0.05', '(5,0.05)', rbf_kernel(X, 5, 0.05, 0.1))
datagen.list[[3]] <- dgobj('rbf10-0.05', '(10,0.05)', rbf_kernel(X, 10, 0.05, 0.1))
png(paste0(save.slug, 'rbf-005-fz01.png'), height=6, width=6, units='in', res = 100)
plot_nice(2, datagen.list, c(-6,6))
dev.off()
 

n.plot <- 4
N <- 100
X <- seq(-1, 1, length.out=N)
K <- rbf_kernel(X, 1, 1, 1e-8)

Y.plot <- list() 
for (i in 1:n.plot) {
  err <- t(mvtnorm::rmvnorm(1, rep(0, N), K))
  mu <- rep(0, N)
  Y.plot[[i]] <- Y <- mu + err
}

x.restrict <- 1:25

png('./package/gp_kernlearn/code/output/acv/data_gen_eg_restricted.png', height=5, width=5, units='in', res=100)
plot(NA, NA, xlim=c(-1,1), ylim=c(-2,2), xlab='x', ylab='')
for (i in 1:n.plot) {
  lines(X[x.restrict], Y.plot[[i]][x.restrict])
}
dev.off()

png('./package/gp_kernlearn/code/output/acv/data_gen_eg.png', height=5, width=5, units='in', res=100)
plot(NA, NA, xlim=c(-1,1), ylim=c(-2,2), xlab='x', ylab='')
for (i in 1:n.plot) {
  lines(X, Y.plot[[i]])
}
dev.off()




rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/kernlearn.R')

mod.l5.rbf <- readRDS('package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_tr100_legd5.RData')
X <- mod.l5.rbf$b.X.train[,2]
Y.obs <- mod.l5.rbf$Y.train[[length(mod.l5.rbf$Y.train)]]

mod.l7.rbf <- readRDS('package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_tr100_legd7.RData')

l5.post <- posterior_test_meanvar_brev(mod.l5.rbf$E.est, mod.l5.rbf$V.est, matrix(X, ncol=1), Y.obs, matrix(X, ncol=1), make_legendre1D_basis_maker(5))
l7.post <- posterior_test_meanvar_brev(mod.l7.rbf$E.est, mod.l7.rbf$V.est, matrix(X, ncol=1), Y.obs, matrix(X, ncol=1), make_legendre1D_basis_maker(7))


png('./package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_fits.png', height=5, width=5, units='in', res=100)
plot(X, Y.obs, pch=16,
     ylim=c(min(Y.obs), max(Y.obs)*1.5), ylab='')
lines(X, l5.post$post.test.mean, lwd=2, col='blue')
#polygon(c(X, rev(X)), c(l5.post$post.test.mean+2*sqrt(diag(l5.post$post.test.var)), rev(l5.post$post.test.mean-2*sqrt(diag(l5.post$post.test.var)))), col='blue', border = NA)
lines(X, l7.post$post.test.mean, lwd=2, col='green')
#polygon(c(X, rev(X)), c(l7.post$post.test.mean+2*sqrt(diag(l7.post$post.test.var)), rev(l7.post$post.test.mean-2*sqrt(diag(l7.post$post.test.var)))), col='green', border = NA)
legend('top', legend=c('Obs', '5th Degree', '7th Degree'), col=c('black','blue','green'), pch=c(16, NA, NA), lty=c(NA, 1, 1), lwd=c(NA, 3,3))
dev.off()




png('./package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_STANfits.png', height=5, width=5, units='in', res=100)
plot(X, Y.obs, pch=16,
     ylim=c(min(Y.obs), max(Y.obs)*1.5), ylab='')
lines(X, mod.l5.rbf$stan.output$fitted.values, lwd=2, col='blue')
lines(X, mod.l7.rbf$stan.output$fitted.values, lwd=2, col='green')
legend('top', legend=c('Obs', '5th Degree', '7th Degree'), col=c('black','blue','green'), pch=c(16, NA, NA), lty=c(NA, 1, 1), lwd=c(NA, 3,3))
dev.off() 

n <- length(X)


png('./package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_Ktrue.png', height=5, width=5, units='in', res=100)
image(mod.l7.rbf$K.true[n:1,], xaxt='n', yaxt='n')
dev.off()
png('./package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_tr100_legd5_Kest.png', height=5, width=5, units='in', res=100)
image(mod.l7.rbf$K.est[n:1,], xaxt='n', yaxt='n')
dev.off()
png('./package/gp_kernlearn/code/output/gof_bs/models/model_Yrbf5-005-fz01_tr100_legd7_Kest.png', height=5, width=5, units='in', res=100)
image(mod.l5.rbf$K.est[n:1,], xaxt='n', yaxt='n')
dev.off()
