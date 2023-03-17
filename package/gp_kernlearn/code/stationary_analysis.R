rm(list=ls())

setwd('~/Git/gp_kernlearn/')
#source('package/gp_kernlearn/code/basis_orthog_poly.R')
#source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/stationary.R')

get_unique <- function(D, K) {
  
}

save_slug <- 'package/gp_kernlearn/code/output/stat/'
dir.create(file.path(save_slug))

model.name <- 'model_Yrbf5-05-fz01_tr100_legd10'
saved.model.path <- paste0('package/gp_kernlearn/code/output/gof_bs/models/', model.name, '.RData')

model <- readRDS(saved.model.path)  
X <- model$b.X.train[,2] # TODO it would be cool if kernlearn just saved X untouched.

#s.x <- s.out.lm$D
#s.y <- s.out.lm$modtrain
K.corr <- model$K.est/max(model$K.est)
s.out.lm <- make_stationary(K.corr, X, dist_euc, modfit_lm)
s.out.p3 <- make_stationary(K.corr, X, dist_euc, modfit_lm3)
s.out.mon3 <- make_stationary(K.corr, X, dist_euc, modfit_monopoly_maker(3))

plot(s.out.lm$D, s.out.lm$modtrain, col='black', cex=0.1, ylab='Correlation', xlab='Distance')
lines(s.out.lm$D[4951:5050], model$K.true[100,]/max(model$K.true), pch=16, col='green', lwd=3)
lines(s.out.lm$D[4951:5050], s.out.lm$K.stationary[100,], col='red', lwd=3)
lines(s.out.p3$D[4951:5050], s.out.p3$K.stationary[100,], col='orange', lwd=3)
lines(s.out.mon3$D[4951:5050], s.out.mon3$K.stationary[100,], col='magenta', lwd=3)
legend('topright', c('Dist(K hat)', 'True', 'Linear', 'Poly 3', 'Mono 3'), 
       col=c('black', 'green', 'red', 'orange', 'magenta'), pch=c(16, rep(NA, 4)), lwd = c(NA, rep(3, 4)))



model$K.true[1:3]
model$K.est[1:3]
s.out.mon3$K.stationary[1:3]

image(model$K.true)
image(model$K.est)
image(s.out.lm$K.stationary)
image(s.out.p3$K.stationary)
image(s.out.mon3$K.stationary)
