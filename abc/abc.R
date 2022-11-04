# Copied from my_abc.R 2023/11/03


library(foreach)
library(doParallel)
library(dplyr)

abc_knn_fixedrt <- function(M, runtime, param_prior, data_gen, y_0, distance, n_cores=1, k=1, packages_load=c()) {
  # Two layer hiearchical model: theta ~ prior; y|theta ~ data_gen
  #
  # fixed number of samples returned, fixed runtime, KNN for acceptance
  # for theta p-dim, d(s_0, s') k dim, this returns 
  # M x (p+k) matrix of simulated theta and accompanying distances.  
  # runtime: minutes
  
  
  if (is.null(dim(y_0))) {
    # 1D response 
    n <- length(y_0)
  } else {
    # vector response
    n <- dim(y_0)[1]
  }
  
  if (k != 1) {
    errorCondition('have not implemented knn for higher dim distance metrics')
  }
  
  # theta in Rp
  p <- length(unlist(param_prior()))
  
  runtime <- 60*runtime
  
  cl<-makeCluster(n_cores)
  registerDoParallel(cl)

  abc.output <- foreach (c = 1:n_cores, .packages=packages_load) %dopar% {
    start <- proc.time()[3]
    all.sim <- as.data.frame(matrix(NA, nrow=0, ncol=(p+k)))
    names(all.sim) <- 1:(p+k)
    
    while (proc.time()[3] - start < runtime) {
      theta <- param_prior()
      y <- data_gen(n, theta)
      d <- distance(y, y_0)
      
      new <-  as.data.frame(matrix(c(unlist(theta), d), nrow=1, dimnames=list(NA, 1:(p+k))))
      all.sim <- dplyr::bind_rows(all.sim, new) # dplyr not crazy slow
    }
    
    all.sim
  }
  stopCluster(cl)
  
  abc.output <- bind_rows(abc.output)
  abc.output <- abc.output[order(abc.output[,p+k]),] # final col is distance metric, TODO: brittle 
  colnames(abc.output) <- c(names(unlist(param_prior())), 'dist')
  
  returned <- list(n_gen=nrow(abc.output),
                   results=NULL)

  if (nrow(abc.output) > M) {
    returned$results <- abc.output[1:M,]
  } else {
    cat('Warning: less than M samples generated in alloted time')
    returned$results <- abc.output # columms: theta_1, theta2,...,theta_p, dist1, ..distk
  }

  return(returned)
}




abc_thresh_fixedM <- function(M, eps, param_prior, data_gen, y_0, distance, n_cores=1, k=1, packages_load=c()) {
  # Two layer hiearchical model: theta ~ prior; y|theta ~ data_gen
  #
  # fixed number of samples M returned, variable runtime, hard threshold for acceptance
  # for theta p-dim, d(s_0, s') k dim, returns 
  # n x (p+k) matrix of simulated theta and accompanying distances.  
  
  if (is.null(dim(y_0))) {
    # 1D response 
    n <- length(y_0)
  } else {
    # vector response
    n <- dim(y_0)[1]
  }
  
  if (k != 1) {
    errorCondition('have not implemented knn for higher dim distance metrics')
  }
  
  # theta in Rp
  p <- length(unlist(param_prior()))
  
  
  start <- proc.time()[3]
  cl<-makeCluster(n_cores)
  registerDoParallel(cl)
  
  batch.sizes <- rep(floor(M/n_cores), n_cores)
  if (n_cores > 1 ) {
    batch.sizes[n_cores] <- M - floor(M/n_cores)*(n_cores-1)  
  }
  
  abc.output <- foreach (c = 1:n_cores, .packages=packages_load) %dopar% {
    out <- as.data.frame(matrix(NA, nrow=batch.sizes[c], ncol=(p+k)))
    names(out) <- 1:(p+k)
    
    accepted <- 0
    while (accepted < batch.sizes[c]) {
      theta <- param_prior()
      y <- data_gen(n, theta)
      d <- distance(y, y_0)
      
      if (d <= eps) {
        accepted <- accepted + 1
        out[accepted, ] <- c(unlist(theta), d)
      }
      #new <-  as.data.frame(matrix(c(unlist(theta), d), nrow=1, dimnames=list(NA, 1:(p+k))))
      #all.sim <- dplyr::bind_rows(all.sim, new) # dplyr not crazy slow
      
    }
    
    out
  }
  stopCluster(cl)
  
  abc.output <- bind_rows(abc.output)
  
  if (M != nrow(abc.output)) {
    stop('something is wrong; returning ', nrow(abc.output), ' rows.\n')
  }
  
  #abc.output <- abc.output[order(abc.output[,p+k]),] # final col is distance metric, TODO: brittle 
  colnames(abc.output) <- c(names(unlist(param_prior())), 'dist')
  
  returned <- list(runtime=proc.time()[3] - start,
                   results=abc.output)
  
  return(returned)
  
  # out <- as.data.frame(matrix(NA, nrow=0, ncol=(p+k)))
  # names(out) <- 1:(p+k)
  # accepted <- 0 
  # 
  # while (accepted < M) {
  #   theta <- param_prior()
  #   y <- data_gen(n, theta)
  #   d <- distance(y, y_0)
  #   
  #   if (d <= eps) {
  #     new <-  as.data.frame(matrix(c(unlist(theta), d), nrow=1, dimnames=list(NA, 1:(p+k))))
  #     out <- bind_rows(out, new) # dplyr not crazy slow
  #     accepted = 1+accepted
  #   }
  # }
  
  #return(out)
}



abc_thresh_fixedrt <- function(runtime, eps, param_prior, data_gen, y_0, distance, n_cores=1, k=1, packages_load=c()) {
  # Two layer hiearchical model: theta ~ prior; y|theta ~ data_gen
  #
  # fixed computation budget, random number of samples returned, hard threshold for acceptance
  # for theta p-dim, d(s_0, s') k dim, returns 
  # n x (p+k) matrix of simulated theta and accompanying
  # distances.  

  if (is.null(dim(y_0))) {
    # 1D response 
    n <- length(y_0)
  } else {
    # vector response
    n <- dim(y_0)[1]
  }
  
  if (k != 1) {
    errorCondition('have not implemented knn for higher dim distance metrics')
  }
  
  # theta in Rp
  p <- length(unlist(param_prior()))
  
  runtime <- 60*runtime
  
  cl<-makeCluster(n_cores)
  registerDoParallel(cl)
  
  abc.output <- foreach (c = 1:n_cores, .packages=packages_load) %dopar% {
    start <- proc.time()[3]
    out <- as.data.frame(matrix(NA, nrow=0, ncol=(p+k)))
    names(out) <- 1:(p+k)
    
    while (proc.time()[3] - start < runtime) {
      theta <- param_prior()
      y <- data_gen(n, theta)
      d <- distance(y, y_0)
      
      if (d <= eps) {
        new <-  as.data.frame(matrix(c(unlist(theta), d), nrow=1, dimnames=list(NA, 1:(p+k))))
        out <- dplyr::bind_rows(out, new) # dplyr not crazy slow
      }
    }
    
    out
  }
  stopCluster(cl)
  
  abc.output <- bind_rows(abc.output)
  
  return(list(M=nrow(abc.output),
              results=abc.output))
  
  # out <- as.data.frame(matrix(NA, nrow=0, ncol=(p+k)))
  # names(out) <- 1:(p+k)
  # start <- proc.time()[3]
  # 
  # while (proc.time()[3] - start < runtime) {
  #   theta <- param_prior()
  #   y <- data_gen(n, theta)
  #   d <- distance(y, y_0)
  #   
  #   if (d <= eps) {
  #     new <-  as.data.frame(matrix(c(unlist(theta), d), nrow=1, dimnames=list(NA, 1:(p+k))))
  #     out <- bind_rows(out, new) # dplyr not crazy slow
  #   }
  # }
  # 
  # return(out)
}




abc_regression_forward <- function(runtime, eps, param_prior, data_gen, y_0, distance, p, k) {
  # Three layer hierarchical model: 
}





