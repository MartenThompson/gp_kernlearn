rm(list=ls())

setwd('~/Git/gp_kernlearn/')
#source('package/gp_kernlearn/code/basis_orthog_poly.R')
#source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/stationary.R')

#save_slug <- 'package/gp_kernlearn/code/output/stat/'
#dir.create(file.path(save_slug))


model.name <- 'model_Yrbf5-05-fz01_tr100_legd10'
saved.model.path <- paste0('package/gp_kernlearn/code/output/gof_bs/models/', model.name, '.RData')

model <- readRDS(saved.model.path)  
X <- model$b.X.train[,2] # TODO it would be cool if kernlearn just saved X untouched.

#s.x <- s.out.lm$D
#s.y <- s.out.lm$modtrain
s.out.lm <- make_stationary(model$K.est, X, dist_euc, modfit_lm)
s.out.p3 <- make_stationary(model$K.est, X, dist_euc, modfit_lm3)
s.out.mon3 <- make_stationary(model$K.est, X, dist_euc, modfit_monopoly_maker(3))

plot(s.out.lm$D, s.out.lm$modtrain, col='grey')
points(s.out.lm$D, model$K.true[upper.tri(model$K.true, diag=TRUE)], pch=16, col='green')
points(s.out.lm$D, s.out.lm$modpred, pch=16, col='red')
points(s.out.p3$D, s.out.p3$modpred, pch=16, col='orange')
points(s.out.mon3$D, s.out.mon3$modpred, pch=16, col='magenta')
legend('topright', c('Dist(K hat)', 'Linear', 'Poly 3', 'Mono 3'), 
       col=c('grey', 'green', 'red', 'orange', 'magenta'), pch = 16)


model$K.true[1:3]
model$K.est[1:3]
s.out.mon3$K.stationary[1:3]

image(model$K.true)
image(model$K.est)
image(s.out.lm$K.stationary)
image(s.out.p3$K.stationary)
image(s.out.mon3$K.stationary)
