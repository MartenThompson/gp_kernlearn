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
 


