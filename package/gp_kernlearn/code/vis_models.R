rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
#saved.model.path <- paste0('package/gp_kernlearn/code/output/gof_bs/models/', model.name, '.RData')


dgobj <- function(name, legend, K) {
  return(list(name=name, 
              legend=legend,
              K=K))
}

plot_nice <- function() {
  
}

N.x <- 100         # number of x locations
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)

datagen.list <- list()
datagen.list[[1]] <- dgobj('lin1-1-fz01', lin_kernel(X, 1, 1, 0.1))
datagen.list[[2]] <- dgobj('quad1-1-1-fz01', quad_kernel(X,1, 1, 1, 0.1))
datagen.list[[3]] <- dgobj('cubic', cubic_kernel(X, 1,1,1,1,0.1))
  
plot.colors <- viridis::magma(5,0.75)

plot(NA, NA, xlim=c(-1,1), ylim=c(-3,3))
n.viz <- 1
for (i in 1:length(datagen.list)) {
  for (r in 1:n.viz) {
    Y <- mvtnorm::rmvnorm(1, sigma=datagen.list[[i]]$K)
    lines(X, Y, lwd=3, col=plot.colors[i])  
  }
}
legend('topright', legend=lapply(datagen.list, function(l){l$name}), col = plot.colors, lty=1, lwd=3)
