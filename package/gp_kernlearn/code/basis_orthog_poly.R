#### TITLE ####


##############################
#### KERN LEARN FUNCTIONS ####
##############################



#library(orthopolynom)
#temp <- legendre.polynomials(3, normalized = TRUE)


#leg2 <- 
#leg3 <- 

legendre_polynomials <- c( # 0 is identity
  function(x) { # 1
    return(x)
  }, 
  function(x) { # 2
    return(0.5*(3*x^2 -1))
  }, 
  function(x) { # 3
    return(0.5*(5*x^3 - 3*x))
  },
  function(x) { # 4
    return((1/8)*(35*x^4-30*x^2 + 3))
  },
  function(x) { # 5
    return((1/8)*(63*x^5 - 70*x^3 +15*x))
  })


# basis_maker()s should take one argument, X, and have everything else configured.

make_legendre1D_basis_maker <- function(degree) {
  if (5 < degree) {
    stop('Only implemented up to degree 5.')
  }
  
  basis_maker <- function(X) {
    p <- dim(X)[2]
    if (1 != p) {
      stop('X should be Nx1; is Nx', p)
    }
    
    n <- nrow(X)
    x <- X[,1]
    X.basis <- matrix(NA, nrow=n, ncol=degree+1)  
    # Legendre 0
    X.basis[,1] <- rep(1,n)
    
    for (deg in 1:degree) {
      f <- legendre_polynomials[[deg]]
      #print(f)
      X.basis[,deg+1] <- f(x)
    }
    
    return(X.basis)
  }
  
  return(basis_maker)
}


testing <- FALSE
if(testing) {
  f <- legendre_polynomials[[1]]
  x <- seq(-1,1,length.out=10)
  cat(all(x == f(x)))
  
  basis_maker <- make_legendre1D_basis_maker(4)
  X.b <- basis_maker(matrix(x, ncol=1))
  cat(dim(X.b)[1] ==10)
  cat(dim(X.b)[2] ==5)
}


make_legendre_design_matrix_2D <- function(degree, X_normalized) {
  # X_normalized should have two columns in [-1,1]
  p <- 6
  n <- nrow(X_normalized)
  temp <- X_normalized[,1]
  conc <- X_normalized[,2]
  X_dm <- matrix(NA, nrow=n, ncol=p)  
  # Legendre 0
  X_dm[,1] <- rep(1,n)
  # Legendre 1
  X_dm[,2] <- leg1(temp)
  X_dm[,3] <- leg1(conc)
  # Legendre 2
  X_dm[,4] <- leg1(temp)*leg1(conc)
  X_dm[,5] <- leg2(temp)
  X_dm[,6] <- leg2(conc)
  
  return(X_dm)
}











#############################
## stale code ##
# temp.test <- seq(min(X.train[,1]), max(X.train[,1]), length.out=n.test)
# conc.test <- seq(min(X.train[,2]), max(X.train[,2]),length.out=n.test)
# X.test <- make_legendre_design_matrix(NA, expand.grid(temp.test, conc.test, KEEP.OUT.ATTRS = FALSE))
# 
# 
# pred.data <- data.frame(
#   temp=rep(rep(temp.test, n.test), 3),
#   conc=rep(rep(conc.test, each=n.test), 3),
#   gibb=c(test.post.mean, test.post.mean+2*sqrt(diag(test.post.var)), test.post.mean-2*sqrt(diag(test.post.var))),
#   obs_pred=c(rep('predicted', length(test.post.mean)), rep('2 SE', 2*length(test.post.mean)))
# )
# # or 
# pred.data <- data.frame(
#   temp = rep(temp.test, n.test),
#   conc = rep(conc.test, each=n.test),
#   gibb = test.post.mean,
#   obs_pred = rep('predicted', length(temp.test))
# )
# 
# plot.data <- rbind(pred.data,
#                    data.frame(temp=dat$X[,1],
#                         conc=dat$X[,2],
#                         gibb=dat$Y,
#                         obs_pred=rep('obs', n)))
# 
# 
# plot_ly(x=plot.data$temp, y=plot.data$conc, z=plot.data$gibb, type="scatter3d", mode="markers",color=plot.data$obs_pred)

