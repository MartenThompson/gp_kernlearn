rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/norms.R')
source('package/gp_kernlearn/code/kernlearn.R')

save_slug <- 'package/gp_kernlearn/code/output/gof_bs/obs_dists/'
dir.create(file.path(save_slug))


gen_data <- function(n.train.datasets, N.x, Kern) {
  Y.train <- list()
  for (r in 1:n.train.datasets) {
    err <- t(mvtnorm::rmvnorm(1, rep(0, N.x), Kern))
    mu <- rep(0, N.x)
    Y <- mu + err
    Y.train[[r]] <- Y
  }
  obs.mat <- matrix(unlist(Y.train), ncol=N.x, byrow = TRUE)
  return(obs.mat)  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Training Data Generation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2023)

N.x <- 100         # number of x locations
n.train.datasets <- 100
X <- matrix(seq(-1,1,length.out=N.x), nrow=N.x, ncol=1)

K.lin <- lin_kernel(X, 1, 1, 0.5)
K.quad <- quad_kernel(X, 1, 1, 1, 0.5)
K.cubic <- cubic_kernel(X, 1, 1, 1, 1, 0.5)

e.lin <- eigen(K.lin, symmetric = T)
e.quad <- eigen(K.quad, symmetric = T)
e.cubic <- eigen(K.cubic, symmetric = T)

e.lin$values[1:5]
e.quad$values[1:5]
e.cubic$values[1:5]

plot(log(ts(e.lin$values[1:7])))
lines(log(ts(e.quad$values[1:7])))
lines(log(ts(e.cubic$values[1:7])))

obs.lin <- gen_data(n.train.datasets, N.x, K.lin)
obs.quad <- gen_data(n.train.datasets, N.x, K.quad)
obs.cubic <- gen_data(n.train.datasets, N.x, K.cubic)

e.obs.lin <- eigen(cov(obs.lin), symmetric = T)
e.obs.quad <- eigen(cov(obs.quad), symmetric = T)
e.obs.cubic <- eigen(cov(obs.cubic), symmetric = T)

e.obs.lin$values[1:5]
e.obs.quad$values[1:5]
e.obs.cubic$values[1:5]

plot(log(ts(e.obs.lin$values[1:7])))
lines(log(ts(e.obs.quad$values[1:7])))
lines(log(ts(e.obs.cubic$values[1:7])))
# obvious elbow for all 3 when n.x = 100, n.train.ds = 1e3
# and when n.x = 1d3, n.train.ds = 100 (slow tho)
#image(K.quad)
#image(K.cubic)


## ##
K.rbf1 <- rbf_kernel(X, 5, 5, 0.1)
K.rbf2 <- rbf_kernel(X, 5, .5, 0.1)
K.rbf3 <- rbf_kernel(X, 5, .05, 0.1)

e.rbf1 <- eigen(K.rbf1, symmetric = T)
e.rbf2 <- eigen(K.rbf2, symmetric = T)
e.rbf3 <- eigen(K.rbf3, symmetric = T)


plot(log(ts(e.rbf1$values[1:30])))
lines(log(ts(e.rbf2$values[1:30])))
lines(log(ts(e.rbf3$values[1:30])))

obs.rbf1 <- gen_data(n.train.datasets, N.x, K.rbf1)
obs.rbf2 <- gen_data(n.train.datasets, N.x, K.rbf2)
obs.rbf3 <- gen_data(n.train.datasets, N.x, K.rbf3)

e.obs.rbf1 <- eigen(cov(obs.rbf1), symmetric = T)
e.obs.rbf2 <- eigen(cov(obs.rbf2), symmetric = T)
e.obs.rbf3 <- eigen(cov(obs.rbf3), symmetric = T)

plot(log(ts(e.obs.rbf1$values[1:30])))
lines(log(ts(e.obs.rbf2$values[1:30])))
lines(log(ts(e.obs.rbf3$values[1:30])))


# Moral of the story (only supported by one run mind you): eigen values can 
# differentiate between lin, quad, cubic, and different flavors of RBF. This is 
# true for both the data-generating K's themselves and cov(obs.mat).
# So, can we use this to develop a good summary stat to bootstrap (thus keeping
# our overall setup)? 
# Or, should we think of something entirely different?

## central differences (proxy for 2nd deriv) ##
# use location of max(central_difference) as sum stat?
plot(ts(central_difference(log(e.lin$values[1:7]))), col='red')
lines(ts(central_difference(log(e.quad$values[1:7]))), col='blue')
lines(ts(central_difference(log(e.cubic$values[1:7]))))

plot(ts(central_difference(log(e.obs.lin$values[1:7]))), col='red')
lines(ts(central_difference(log(e.obs.quad$values[1:7]))), col='blue')
lines(ts(central_difference(log(e.obs.cubic$values[1:7]))))


plot(ts(central_difference(log(e.obs.rbf1$values[1:30]))), col='red')
lines(ts(central_difference(log(e.obs.rbf2$values[1:30]))), col='blue')
lines(ts(central_difference(log(e.obs.rbf3$values[1:30]))))


# https://stackoverflow.com/questions/4471993/compute-the-elbow-for-a-curve-automatically-and-mathematically
# https://raghavan.usc.edu//papers/kneedle-simplex11.pdf
