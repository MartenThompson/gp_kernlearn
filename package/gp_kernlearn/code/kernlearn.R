kern <- function(B,V) {
  n <- nrow(B)
  K <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      K[i,j] <- B[i,]%*%V%*%t(t(B[j,]))
    }
  }
  
  K[lower.tri(K)] <- t(K)[lower.tri(K)]
  diag(K) <- diag(K) + 1e-4
  return(K)
}


#######################
#### ABC FUNCTIONS ####
#######################

# TODO make a maker w/ p considered
param_prior_maker <- function(basis.dim, sd=2) {
  param_prior <- function() {
    return(list(beta=rnorm(basis.dim,0,sd)))
  }
  
  return(param_prior)
}

data_gen_maker <- function(X_dm) {
  data_gen <- function(n_depricated, params) {
    return(X_dm%*%params$beta)
  }
  return(data_gen)
}

distance <- function(y, y_0) {
  return(sum((y-y_0)^2))
}


######################################
#### ABC POSTERIOR APPROXIMATIONS ####
######################################

#files <- list.files('ms_eg/jan2022_bcc_files_trimmed/')

#for (f in files) {
#  dat <- read_normalize_unit_square(f)
#  y_0 <- dat$Y
#  n <- length(y_0)
#  X_dm <- make_legendre_design_matrix(NA, dat$X) # this is really b(x)'s
#  data_gen <- data_gen_maker(X_dm)
#  garbage <- data_gen(NA, param_prior()) # just need to run this method once o/w lazy init bites us.
#  
#  M = 100
#  runtime <- 10
#  abc_output <- abc_knn_fixedrt(M, runtime, param_prior, data_gen, y_0, distance, n_cores=6, k=1, packages_load=c())
#  n_gen <- abc_output$n_gen
#  dataset_name <- strsplit(f, '.', fixed=TRUE)[[1]][1]
#  abc_output_name <- paste0('ms_eg/output/abc_output/', dataset_name,
#                            '_', runtime, 'min6_l2_M100n', n_gen, '.csv')
#  write.csv(abc_output$results, abc_output_name, row.names = FALSE)
#}




# An intelligent thing to do would be incorporate ABC uncertainty in a noise kernel.

#par(mfrow=c(2,3))
#for (i in 1:6) {
#  plot(density(abc_output$results[,i]), main=paste0('Beta ', i))
#  lines(c(E.leg[i], E.leg[i]), c(0,1), col='red')
#}
#par(mfrow=c(1,1))

# for (f in files) {
#   dat <- read_normalize_unit_square(f)
#   temp.midpoint <- max(dat$X[,1]) - (max(dat$X[,1]) - min(dat$X[,1]))/2
#   train_on_lowT <- dat$X[,1] < temp.midpoint
#   
#   y_0 <- dat$Y[train_on_lowT]
#   n <- length(y_0)
#   X_dm <- make_legendre_design_matrix(NA, dat$X[train_on_lowT,]) # this is really b(x)'s
#   data_gen <- data_gen_maker(X_dm)
#   garbage <- data_gen(NA, param_prior()) # just need to run this method once o/w lazy init bites us.
#   
#   M = 100
#   runtime <- 10
#   abc_output <- abc_knn_fixedrt(M, runtime, param_prior, data_gen, y_0, distance, n_cores=6, k=1, packages_load=c())
#   n_gen <- abc_output$n_gen
#   dataset_name <- strsplit(f, '.', fixed=TRUE)[[1]][1]
#   abc_output_name <- paste0('ms_eg/output/abc_output/', dataset_name,
#                             'trainlowT_', runtime, 'min6_l2_M100n', n_gen, '.csv')
#   write.csv(abc_output$results, abc_output_name, row.names = FALSE)
# }




##########################
### MAKING THE KERNEL ####
##########################

make_kern <- function (X.train, X.test, V.basis, basis_maker) {
  # n.test <- 30 # it will be this squared fyi  
  n <- nrow(X.train)
  X.trte <- rbind(X.train, X.test)
  X.basis.trte <- basis_maker(X.trte)
  K.trte <- kern(X.basis.trte, V.basis)
  #image(K.trte)
  
  
  K.train <- K.trte[1:n, 1:n]
  K.trstar <- K.trte[1:n,(n+1):ncol(K.trte)]
  K.startr <- K.trte[(n+1):ncol(K.trte), 1:n]
  K.star <- K.trte[(n+1):ncol(K.trte), (n+1):ncol(K.trte)]
  
  return(list(
    K.train=K.train,
    K.traininv = solve(K.train),
    K.startr=K.startr,
    K.star=K.star
  ))
}


posterior_test_meanvar <- function(abc.samples, X.train, Y.train, X.test, basis_maker) {
  V <- cov(as.matrix(abc.samples)) 
  E <- apply(abc.samples, 2, mean) 
  
  kern.pieces <- make_kern(X.train, X.test, V, basis_maker)
  K.traininv <- kern.pieces$K.traininv
  K.startr <- kern.pieces$K.startr
  K.star <- kern.pieces$K.star
  
  
  X.test.basis <- basis_maker(X.test) 
  X.train.basis <- basis_maker(X.train)
  
  post.test.mean <- X.test.basis%*%E + K.startr%*%K.traininv%*%(Y.train - X.train.basis%*%E)
  post.test.var <- K.star - K.startr%*%K.traininv%*%t(K.startr)
  
  return(list(
    post.test.mean = post.test.mean,
    post.test.var = post.test.var,
    kern.pieces = kern.pieces
  ))
}

posterior_test_meanvar_br <- function(rstan_output, X.train, Y.train, X.test, basis_maker) {
  E <- unname(rstan_output$coefficients)
  V <- unname(rstan_output$covmat)
  
  kern.pieces <- make_kern(X.train, X.test, V, basis_maker)
  K.traininv <- kern.pieces$K.traininv
  K.startr <- kern.pieces$K.startr
  K.star <- kern.pieces$K.star
  
  
  X.test.basis <- basis_maker(X.test) 
  X.train.basis <- basis_maker(X.train)
  
  post.test.mean <- X.test.basis%*%E + K.startr%*%K.traininv%*%(Y.train - X.train.basis%*%E)
  post.test.var <- K.star - K.startr%*%K.traininv%*%t(K.startr)
  
  return(list(
    post.test.mean = post.test.mean,
    post.test.var = post.test.var,
    kern.pieces = kern.pieces
  ))
}

posterior_test_meanvar_brev <- function(E, V, X.train, Y.train, X.test, basis_maker) {
  kern.pieces <- make_kern(X.train, X.test, V, basis_maker)
  K.traininv <- kern.pieces$K.traininv
  K.startr <- kern.pieces$K.startr
  K.star <- kern.pieces$K.star
  
  
  X.test.basis <- basis_maker(X.test) 
  X.train.basis <- basis_maker(X.train)
  
  post.test.mean <- X.test.basis%*%E + K.startr%*%K.traininv%*%(Y.train - X.train.basis%*%E)
  post.test.var <- K.star - K.startr%*%K.traininv%*%t(K.startr)
  
  return(list(
    post.test.mean = post.test.mean,
    post.test.var = post.test.var,
    kern.pieces = kern.pieces
  ))
}

####################
##  TOY KERNELS  ##
####################
# assumes 0 centered X
lin_kernel <- function(X, sigma.0, sigma.1, fuzz=0.1) {
  K <- sigma.0 + sigma.1*X%*%t(X)
  diag(K) <- diag(K) + fuzz
  return(K)
}

quad_kernel <- function(X, sigma.0, sigma.1, sigma.2, fuzz=0.1) {
  #K <- sigma.0 + sigma.1*X%*%t(X) + sigma.2*(X%*%t(X))%*%(X%*%t(X))
  #K <- (sigma.0 + sigma.1*X%*%t(X))%*%(sigma.0 + sigma.2*X%*%t(X))
  n <- dim(X)[1]
  K <- matrix(NA, n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- sigma.0 + sigma.1*X[i,]%*%t(X[j,]) + sigma.2*X[i,]%*%t(X[j,])%*%X[i,]%*%t(X[j,])
    }
  }
  diag(K) <- diag(K) + fuzz
  return(K)
}

cubic_kernel <- function(X, sigma.0, sigma.1, sigma.2, sigma.3, fuzz=0.1) {
  n <- dim(X)[1]
  K <- matrix(NA, n,n)
  for (i in 1:n) {
    for (j in 1:n) {
      K[i,j] <- sigma.0 + sigma.1*X[i,]%*%t(X[j,]) + 
        sigma.2*X[i,]%*%t(X[j,])%*%X[i,]%*%t(X[j,]) +
        sigma.3*X[i,]%*%t(X[j,])%*%X[i,]%*%t(X[j,])%*%X[i,]%*%t(X[j,])
    }
  }
  diag(K) <- diag(K) + fuzz
  return(K)
}

rbf_kernel <- function(X, alpha.0, alpha.1) {
  K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
  return(K)
}

####################
## USING ALL THIS ##
####################
library(rstanarm)

fit_kernlearn <- function(b.X, basis_maker, Y.list) {
  degree <- dim(b.X)[2]-1
  beta.manysamples.pri <- NA
  beta.manysamples.alt <- NA
  n.train.data <- length(Y.list)
  
  for (r in 1:n.train.data) {
    Y <- Y.list[[r]]
    
    silent <- capture.output(br.out <- stan_glm(Y ~ b.X[,2:(degree+1)], family=gaussian()))
    beta.samples <- matrix(unlist(br.out$stanfit@sim$samples[[1]][[1]]), ncol=1)
    for (i in 2:(degree+1)) {
      new <- matrix(unlist(br.out$stanfit@sim$samples[[1]][[i]]), ncol=1)
      beta.samples <- cbind(beta.samples, new)
    }
    
    if (1==r) {
      beta.manysamples <- beta.samples
    } else {
      beta.manysamples <- rbind(beta.manysamples, beta.samples)
    }
  }
  
  E <- apply(beta.manysamples, 2, mean)
  V <- cov(beta.manysamples)
  
  X.test <- matrix(seq(-5,5,length.out=20), ncol=1)
  post.ests <- posterior_test_meanvar_brev(E, V, X, Y.list[[1]], X.test, basis_maker)
  K.hat <- post.ests$kern.pieces$K.train
  
  # # 1D only!
  # model <- function(x) {
  #   x <- matrix(x, ncol=1)
  #   post.ests <- posterior_test_meanvar_brev(E, V, X, Y.list[[1]], X, basis_maker)
  # return mean like function? unclear what y should be then.
  # }
  
  return(list(
    K.est = K.hat,
    stan.output =br.out,
    Y.train = Y.list,
    b.X.train = b.X,
    basis_maker = basis_maker
  ))
}


make_predictions <-function(f, slug) {
  cat(f, '\n')
  dat <- read_normalize_unit_square(f)
  temp.midpoint <- max(dat$X[,1]) - (max(dat$X[,1]) - min(dat$X[,1]))/2
  train_on_lowT <- dat$X[,1] < temp.midpoint
  
  y.0 <- dat$Y[train_on_lowT]
  X.train <- dat$X[train_on_lowT,]
  y.test <- dat$Y[!train_on_lowT]
  X.test <- dat$X[!train_on_lowT,]
  
  dataset_name <- strsplit(f, '.', fixed=TRUE)[[1]][1]
  abc_output_name <- paste0('^', dataset_name, slug)
  abc_output <- read.csv(paste0('ms_eg/output/abc_output/', 
                                dir('ms_eg/output/abc_output/', pattern=abc_output_name)[1]))
  cat(nrow(X.train), ' ', nrow(X.test), '\n')
  V.leg <- cov(as.matrix(abc_output[,1:6]))
  E.leg <- apply(abc_output[,1:6], 2, mean)
  kern_pieces <- make_kern(X.train, X.test, V.leg)
  post_ests <- posterior_test_meanvar(X.test, X.train, kern_pieces$K.traininv, 
                                      kern_pieces$K.startr, kern_pieces$K.star,
                                      E.leg, y.0)
  
  pred_data <- data.frame(
    temp = X.test[,1],
    conc = X.test[,2],
    y.true = y.test,
    pred = post_ests$post.test.mean,
    pred_se = 2*sqrt(diag(post_ests$post.test.var))
  )
  
  obs_data <- data.frame(
    temp = X.train[,1],
    conc = X.train[,2],
    y.true = y.0,
    pred = rep(NA, length(y.0)),
    pred_se = rep(NA, length(y.0))
  )
  
  return(list(
    preds = pred_data,
    obs = obs_data
  ))
}


# itr <- 13
# output <- make_predictions(files[itr], 'trainlowT_10min6_l2_M100n')
# plot.data <- data.frame(
#   temp = c(output$preds$temp, output$preds$temp, output$obs$temp),
#   conc = c(output$preds$conc, output$preds$conc, output$obs$conc),
#   gibb = c(output$preds$pred, output$preds$y.true, output$obs$y.true),
#   obs_or_pred = c(rep('pred', nrow(output$preds)), rep('obs', nrow(output$preds)+nrow(output$obs)))
# )


#plot_ly(x=plot.data$temp, y=plot.data$conc, z=plot.data$gibb, type="scatter3d", mode="markers",color=plot.data$obs_or_pred)
#paste0(strsplit(files[itr], '.', fixed=TRUE)[[1]][1],'trainlowT_10min6_l2_M100n')
