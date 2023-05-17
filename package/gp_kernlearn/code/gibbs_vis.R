plot3D_piped <- function(X.test, Y.test, Y.pred) {
  plot.data <- data.frame(
    X1=rep(X.test[,1],2),
    X2=rep(X.test[,2],2),
    Z=c(Y.test, Y.pred),
    type=c(rep('Test Values', length(Y.test)), rep('Predicted', length(Y.pred)))
  )
  colors <- c('black', 'green')
  
  p <- plot_ly(plot.data, x=~X1, y=~X2, z=~Z, color=~type, colors=colors, opacity=1) %>%
    add_markers() %>%
    layout(legend = list(x = 0.1, y = .75),
           scene=list(
             xaxis=list(title='Temp.'),
             yaxis=list(title='Conc.'),
             zaxis=list(title='Gibbs Energy')))
  
  #p <- plot_ly(colorscale = c('#FFE1A1', '#683531')) %>%
  #  add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.test),
  #              color=I('black'), name='Observed')%>%#c(rep('Observed', length(Y.test)), rep('Estimated', length(Y.pred)))) %>%
  # add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.pred),
  #              color=I('green'), name='Predicted')%>%
  
  
  return(p)
}

plot3D_many <- function(X.list, Y.list, labels=NA, color.names=NA) {
  p <- plot_ly() 
  
  for (i in 1:length(Y.list)) {
    if (!is.na(labels) & 'Post. Mean' == labels[i]) {
      opac <- 0.5
    } else {
      opac <- 1
    }
    if (is.na(color.names)) {
      p <- add_markers(p, x=X.list[[i]][,1], y=X.list[[i]][,2], z=Y.list[[i]], 
                       opacity=opac)
      p <- layout(p, showlegend=FALSE)
    } else {
      p <- add_markers(p, x=X.list[[i]][,1], y=X.list[[i]][,2], z=Y.list[[i]], 
                       color=labels[i], colors=color.names, opacity=opac)
      p <- layout(p, legend = list(x = 0.1, y = .75))
    }
    
    
  }
  
  p <- layout(p, 
              scene=list(
                xaxis=list(title='Temp.'),
                yaxis=list(title='Conc.'),
                zaxis=list(title='Gibbs Energy')))
  
  #p <- plot_ly(colorscale = c('#FFE1A1', '#683531')) %>%
  #  add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.test),
  #              color=I('black'), name='Observed')%>%#c(rep('Observed', length(Y.test)), rep('Estimated', length(Y.pred)))) %>%
  # add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.pred),
  #              color=I('green'), name='Predicted')%>%
  
  
  return(p)
}


plot3D_halfpred <- function(X.test, Y.test, Y.pred, X.train, Y.train) {
  plot.data <- data.frame(
    X1=c(rep(X.test[,1],2), rep(X.train[,1],1)),
    X2=c(rep(X.test[,2],2), rep(X.train[,2], 1)),
    Z=c(Y.test, Y.pred, Y.train),
    type=c(rep('Test Values', length(Y.test)), 
           rep('Predicted', length(Y.pred)),
           rep('Train Values', length(Y.train)))
  )
  colors <- c('black', 'green', 'gray')
  
  p <- plot_ly(plot.data, x=~X1, y=~X2, z=~Z, color=~type, colors=colors, opacity=1) %>%
    add_markers() %>%
    layout(legend = list(x = 0.1, y = .75),
           scene=list(
             xaxis=list(title='Temp.'),
             yaxis=list(title='Conc.'),
             zaxis=list(title='Gibbs Energy')))
  
  #p <- plot_ly(colorscale = c('#FFE1A1', '#683531')) %>%
  #  add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.test),
  #              color=I('black'), name='Observed')%>%#c(rep('Observed', length(Y.test)), rep('Estimated', length(Y.pred)))) %>%
  # add_markers(x=rep(X.test[,1],1), y=rep(X.test[,2],1), z=c(Y.pred),
  #              color=I('green'), name='Predicted')%>%
  
  
  return(p)
}


one_time_data_clean <- function() {
  stop('Already done; saving for record-keeping')
  dat.path <- './package/gp_kernlearn/data/jan2022_bcc_files/'
  save.path <- './package/gp_kernlearn/code/output/gibbs/'
  
  bcc.files <- list.files(dat.path, pattern = '*.CSV')
  n.total <- length(bcc.files) 
  
  y.min <- NA
  y.max <- NA
  temp.min <- NA
  temp.max <- NA
  for (i in 1:n.total) {
    cat('Reading:', i, ' ')
    fn <- bcc.files[i]
    dat.in <- read.csv(paste0(dat.path, fn), sep=' ')
    Y <- dat.in$G.J.mol.
    x.temperature <- dat.in$T.K.
    #x.conc.a <- dat$X.A.
    y.min <- min(y.min, Y, na.rm=TRUE)
    y.max <- max(y.max, Y, na.rm=TRUE)
    temp.min <- min(temp.min, x.temperature, na.rm=TRUE)
    temp.max <- max(temp.max, x.temperature, na.rm=TRUE)
  }
  
  
  Y.scaled <- list()
  Y.unscaled <- list()
  X.scaled <- list()
  X.unscaled <- list()
  for (i in 1:n.total) {
    cat('Scaling:', i, ' ')
    fn <- bcc.files[i]
    dat.in <- read.csv(paste0(dat.path, fn), sep=' ')
    Y.unscaled[[i]] <- dat.in$G.J.mol.
    X.unscaled[[i]] <- cbind(dat.in$T.K., dat.in$X.A.)
    
    Y <- (dat.in$G.J.mol. + (y.max-y.min)/2)/(5*(y.max-y.min))
    x.temperature <- (dat.in$T.K. + (temp.max-temp.min)/2)/(temp.max - temp.min)
    x.conc.a <- dat.in$X.A.
    X <- cbind(x.temperature, x.conc.a)
    
    Y.scaled[[i]] <- Y
    X.scaled[[i]] <- X
  }
  
  write.csv(data.frame(
    ymin = y.min, ymax=y.max, tempmin=temp.min, tempmax=temp.max
  ), paste0(save.path, 'scale_factors.csv'), row.names = FALSE)
  saveRDS(Y.scaled, paste0(save.path, 'Y_scaled.RData'))
  saveRDS(Y.unscaled, paste0(save.path, 'Y_unscaled.RData'))
  saveRDS(X.scaled, paste0(save.path, 'X_scaled.RData'))
  saveRDS(X.unscaled, paste0(save.path, 'X_unscaled.RData'))
}
# TOO EXPENSIVE !
#post.out <- posterior_test_meanvar_brev(E.beta, V.beta, 
#                                        X.train.long, Y.train.long, 
#                                        X.test, basis_maker)
#(p <- plot3D_piped(X.all[[test.itr]], Y.all[[test.itr]], post.out$post.test.mean))
# TOO EXPENSIVE !
#post.out.obs.templow <- posterior_test_meanvar_brev(E.beta, V.beta,
#                                                    rbind(X.train.long, X.test.templow), c(Y.train.long, Y.test.templow),
#                                                    X.test.temphigh, basis_maker)

#plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow$post.test.mean, X.test.templow, Y.test.templow)
