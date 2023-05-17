rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('package/gp_kernlearn/code/kernlearn.R')
source('package/gp_kernlearn/code/basis_orthog_poly.R')
source('package/gp_kernlearn/code/acv_est.R')

require(plotly)


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

plot3D_many <- function(X.list, Y.list) {
  p <- plot_ly() 
  for (i in 1:length(Y.list)) {
    p <- add_markers(p, x=X.list[[i]][,1], y=X.list[[i]][,2], z=Y.list[[i]])
  }
  
  p <- layout(p, legend = list(x = 0.1, y = .75),
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

plot3D_many <- function(X.list, Y.list) {
  p <- plot_ly() 
  for (i in 1:length(Y.list)) {
    p <- add_markers(p, x=X.list[[i]][,1], y=X.list[[i]][,2], z=Y.list[[i]])
  }
  
  p <- layout(p, legend = list(x = 0.1, y = .75),
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

skinny <- function(kernlearn_out) {
  temp <- kernlearn_out
  temp$stan.output <- kernlearn_out$stan.output$stanfit
  return(temp)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data Import ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.path <- './package/gp_kernlearn/code/output/gibbs/'
Y.all <- readRDS(paste0('./package/gp_kernlearn/code/output/gibbs/', 'Y_scaled.Rdata'))
X.all <- readRDS(paste0('./package/gp_kernlearn/code/output/gibbs/', 'X_scaled.Rdata'))

leg.deg <- 2
basis_maker <- make_legendre2D_basis_maker(leg.deg)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### N train = 2 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E dominates, K does not influence enough to make any predictions signif differ
# from post mean.


n.train <- 2 #floor(0.5*n.total)
b.X.train <- list()
X.train.long <- X.all[[1]]
Y.train.long <- Y.all[[1]]
for (i in 1:n.train) {
  cat(i,' ')
  b.X.train[[i]] <- basis_maker(X.all[[i]])
  
  if (1 != i) {
    X.train.long <- rbind(X.train.long, X.all[[i]])
    Y.train.long <- c(Y.train.long, Y.all[[i]])
  }
}

kernlearn.out <- fit_kernlearn_Xunique(b.X.train, basis_maker, Y.all[1:n.train], X.train.long, Y.train.long)

E.beta <- kernlearn.out$E.est
V.beta <- kernlearn.out$V.est

test.itr <- 10
X.test <- X.all[[test.itr]]
Y.test <- Y.all[[test.itr]]
post.out <- posterior_test_meanvar_brev(E.beta, V.beta, 
                                        X.train.long, Y.train.long, 
                                        X.test, basis_maker)
p <- plot3D_piped(X.all[[test.itr]], Y.all[[test.itr]], post.out$post.test.mean)
saveRDS(list(leg.deg=leg.deg,
             n.train=n.train,
             b.X.train=b.X.train,
             kernlearn.out=skinny(kernlearn.out),
             X.train.long=X.train.long,
             Y.train.long=Y.train.long,
             X.test=X.test,
             post.out=post.out,
             plot.preds=p), paste0(save.path, 'analysis_ntr', n.train, '_te', test.itr, '.RData'))

strsplit(bcc.files[test.itr], '.CSV')[[1]]

summary(X.test)
#x.temp.mid <- max(X.test[,1])-(max(X.test[,1])-min(X.test[,1]))/2 # halfway thru domain, not num values
x.temp.q50 <- quantile(X.test[,1], 0.5)
X.test.templow <- X.test[X.test[,1] < x.temp.q50, ]
X.test.temphigh <- X.test[X.test[,1] >= x.temp.q50, ]
Y.test.templow <- Y.test[X.test[,1] < x.temp.q50]
Y.test.temphigh <- Y.test[X.test[,1] >= x.temp.q50]
post.out.obs.templow <- posterior_test_meanvar_brev(E.beta, V.beta,
<<<<<<< HEAD
                                                    rbind(X.train.long, X.test.templow), c(Y.train.long, Y.test.templow),
                                                    X.test.temphigh, basis_maker)
=======
                                                 rbind(X.train.long, X.test.templow), c(Y.train.long, Y.test.templow),
                                                 X.test.temphigh, basis_maker)
>>>>>>> 0739f12bf395c93660be3ce7ca2e1b30ddf1035d

plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow$post.test.mean, X.test.templow, Y.test.templow)


X.test.grid <- as.matrix(expand.grid(seq(min(X.test[,1]),max(X.test[,1]), length.out=30),
<<<<<<< HEAD
                                     seq(min(X.test[,2]),max(X.test[,2]), length.out=30), KEEP.OUT.ATTRS = FALSE))
post.out.obs.templow.grid <- posterior_test_meanvar_brev(E.beta, V.beta,
                                                         X.test.templow, Y.test.templow,
                                                         X.test.grid, basis_maker)
=======
                           seq(min(X.test[,2]),max(X.test[,2]), length.out=30), KEEP.OUT.ATTRS = FALSE))
post.out.obs.templow.grid <- posterior_test_meanvar_brev(E.beta, V.beta,
                                                               X.test.templow, Y.test.templow,
                                                               X.test.grid, basis_maker)
>>>>>>> 0739f12bf395c93660be3ce7ca2e1b30ddf1035d
plot3D_piped(X.test.grid, rep(0.09987,nrow(X.test.grid)), post.out.obs.templow.grid$post.test.mean)






post.out.obs.templow.skinnycond <- posterior_test_meanvar_brev(E.beta, V.beta,
<<<<<<< HEAD
                                                               X.test.templow, Y.test.templow,
                                                               X.test.temphigh, basis_maker)
=======
                                                    X.test.templow, Y.test.templow,
                                                    X.test.temphigh, basis_maker)
>>>>>>> 0739f12bf395c93660be3ce7ca2e1b30ddf1035d

plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow.skinnycond$post.test.mean, X.test.templow, Y.test.templow)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### N train = 10 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.each <- as.numeric(lapply(Y.all, length))

Y.all.ord <- list()
X.all.ord <- list()
rep.itr <- 1
for (i in 1:length(n.each)) {
  idx.small <- which(n.each == sort(n.each)[i])
  if (1 < length(idx.small)) {
    n.equal.length <- length(idx.small)
    idx.small <- idx.small[rep.itr]
    
    if (rep.itr == n.equal.length) {
      rep.itr <- 1
    } else {
      rep.itr <- rep.itr + 1  
    }
    
  } else {
    rep.itr <- 1
  }
  cat(idx.small, ' ')  
  Y.all.ord[[i]] <- Y.all[[idx.small]]
  X.all.ord[[i]] <- X.all[[idx.small]]
}
Y.all.ord <- Y.all.ord[-1]
X.all.ord <- X.all.ord[-1]
as.numeric(lapply(Y.all.ord, length))

n.train <- 50 #floor(0.5*n.total)
b.X.train <- list()
X.train.long <- X.all.ord[[1]]
Y.train.long <- Y.all.ord[[1]]
for (i in 1:n.train) {
  cat(i,' ')
  b.X.train[[i]] <- basis_maker(X.all.ord[[i]])
  
  if (1 != i) {
    X.train.long <- rbind(X.train.long, X.all.ord[[i]])
    Y.train.long <- c(Y.train.long, Y.all.ord[[i]])
  }
}

plot3D_many(X.all.ord[1:n.train], Y.all.ord[1:n.train])

kernlearn.out <- fit_kernlearn_Xunique(b.X.train, basis_maker, Y.all.ord[1:n.train], X.train.long, Y.train.long, 
                                       stan.chains=1, stan.iter=100, EV.only = TRUE)

E.beta <- kernlearn.out$E.est
V.beta <- kernlearn.out$V.est

test.itr <- n.train + 1
X.test <- X.all[[test.itr]]
Y.test <- Y.all[[test.itr]]
# TOO EXPENSIVE !
#post.out <- posterior_test_meanvar_brev(E.beta, V.beta, 
#                                        X.train.long, Y.train.long, 
#                                        X.test, basis_maker)
#(p <- plot3D_piped(X.all[[test.itr]], Y.all[[test.itr]], post.out$post.test.mean))

#saveRDS(list(leg.deg=leg.deg,
#             n.train=n.train,
#             b.X.train=b.X.train,
#             kernlearn.out=skinny(kernlearn.out),
#             X.train.long=X.train.long,
#             Y.train.long=Y.train.long,
#             X.test=X.test,
#             post.out=post.out,
#             plot.preds=p), paste0(save.path, 'analysis_ntr', n.train, '_te', test.itr, '.RData'))

strsplit(bcc.files[test.itr], '.CSV')[[test.itr]]

summary(X.test)
#x.temp.mid <- max(X.test[,1])-(max(X.test[,1])-min(X.test[,1]))/2 # halfway thru domain, not num values
x.temp.q50 <- quantile(X.test[,1], 0.5)
X.test.templow <- X.test[X.test[,1] < x.temp.q50, ]
X.test.temphigh <- X.test[X.test[,1] >= x.temp.q50, ]
Y.test.templow <- Y.test[X.test[,1] < x.temp.q50]
Y.test.temphigh <- Y.test[X.test[,1] >= x.temp.q50]

# TOO EXPENSIVE !
#post.out.obs.templow <- posterior_test_meanvar_brev(E.beta, V.beta,
#                                                    rbind(X.train.long, X.test.templow), c(Y.train.long, Y.test.templow),
#                                                    X.test.temphigh, basis_maker)

#plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow$post.test.mean, X.test.templow, Y.test.templow)


post.out.obs.templow.skinnycond <- posterior_test_meanvar_brev(E.beta, V.beta,
                                                               X.test.templow, Y.test.templow,
                                                               X.test.temphigh, basis_maker)

plot3D_halfpred(X.test.temphigh, Y.test.temphigh, post.out.obs.templow.skinnycond$post.test.mean, X.test.templow, Y.test.templow)
