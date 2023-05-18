#### TITLE ####


##############################
#### KERN LEARN FUNCTIONS ####
##############################



#library(orthopolynom)
#temp <- legendre.polynomials(3, normalized = TRUE)


monomials <- c(
  function(x) {x}, # 1
  function(x) {x^2}, # 2
  function(x) {x^3}, # 3
  function(x) {x^4}, # 4
  function(x) {x^5}, # 5
  function(x) {x^6}, # 6
  function(x) {x^7}, # 7
  function(x) {x^8}, # 8
  function(x) {x^9}, # 9
  function(x) {x^10} # 10
)

legendre_polynomials <- c( # 0 is 1
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
  },
  function(x) { # 6
    return((1/16)*(231*x^6 - 315*x^4 + 105*x^2 - 5))
  },
  function(x) { # 7 
    return((1/16)*(429*x^7 - 693*x^5 + 315*x^3 - 35*x))
  },
  function(x) { # 8
    return((1/128)*(6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35))
  },
  function(x) { # 9
    return((1/128)*(12155*x^9 - 25740*x^7 + 18018*x^5 - 4620*x^3 + 315*x))
  },
  function(x) { # 10
    return((1/256)*(46189*x^10 - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63))
  })

alt_leg <- function(degree, renormalize=FALSE) {
  leg_poly <- function(X) {
    sapply(X, function(x){(1/2^degree)*sum((choose(degree, 0:degree)^2)*(x-1)^(degree-0:degree)*(x+1)^(0:degree))})
  }
  
  if (renormalize) {
    # so integrate(poly_d * poly_d, -1, 1) = 1
    leg_poly_nrm <- function(X) {
      sqrt((2*degree+1)/2)*leg_poly(X)
    }
    return(leg_poly_nrm)
  }
  
  return(leg_poly)
}

# l1 <- alt_leg(1)
# l2 <- alt_leg(2)
# l3 <- alt_leg(3)
# l4 <- alt_leg(4)
# l5 <- alt_leg(5)
# B.alt <- matrix(c(l1(X), l2(X), l3(X), l4(X), l5(X)), ncol=5, byrow=F)
# zapsmall(t(B.alt)%*%B.alt)
# 
# l1.nrm <- alt_leg(1, TRUE)
# l2.nrm <- alt_leg(2, TRUE)
# l3.nrm <- alt_leg(3, TRUE)
# l4.nrm <- alt_leg(4, TRUE)
# l5.nrm <- alt_leg(5, TRUE)
# integrate(function(x){l3.nrm(x)*l3.nrm(x)}, -1 ,1)
# B.alt.nrm <- matrix(c(l1.nrm(X), l2.nrm(X), l3.nrm(X), l4.nrm(X), l5.nrm(X)), ncol=5, byrow=F)
# zapsmall(t(B.alt.nrm)%*%B.alt.nrm)#*mean(diff(X))

# basis_maker()s should take one argument, X, and have everything else configured.

make_legendre1D_basis_maker <- function(degree) {
  if (10 < degree) {
    stop('Only implemented up to degree 10.')
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
    
    if (0 == degree) {
      return(X.basis)
    }
    
    for (deg in 1:degree) {
      f <- legendre_polynomials[[deg]]
      #print(f)
      X.basis[,deg+1] <- f(x)
    }
    
    return(X.basis)
  }
  
  return(basis_maker)
}

# treats intercept (f(x) = 1) as degree zero
make_legendre2D_basis_maker <- function(d1, d2=d1) {
  basis_maker <- function(X) {
    if (2 != ncol(X)) {
      stop('X must be n x 2, actuall has ', ncol(X), ' columns.')
    }
    
    n <- nrow(X)
    x1 <- X[,1]
    x2 <- X[,2]
    
    p <- (d1+1)*(d2+1)
    X.basis <- matrix(NA, nrow=n, ncol=p)
    
    f1.set <- lapply(0:d1, alt_leg)
    f2.set <- lapply(0:d2, alt_leg)
    
    col.fill <- 1
    for (i in 1:(d1+1)) {
      for (j in 1:(d2+1)) {
        f1 <- f1.set[[i]]
        f2 <- f2.set[[j]]
        
        X.basis[,col.fill] <- f1(x1) * f2(x2)
        col.fill <- col.fill+1
      }
    }
    
    return(X.basis)
  }
  
  return(basis_maker)
}

make_monomial1D_basis_maker <- function(degree) {
  if (10 < degree) {
    stop('Only implemented up to degree 10.')
  }
  
  basis_maker <- function(X) {
    p <- dim(X)[2]
    if (1 != p) {
      stop('X should be Nx1; is Nx', p)
    }
    
    n <- nrow(X)
    x <- X[,1]
    X.basis <- matrix(NA, nrow=n, ncol=degree)  
    # Legendre 0
    #X.basis[,1] <- rep(1,n)
    
    if (0 == degree) {
      return(X.basis)
    }
    
    for (deg in 1:degree) {
      f <- monomials[[deg]]
      #print(f)
      X.basis[,deg] <- f(x)
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

