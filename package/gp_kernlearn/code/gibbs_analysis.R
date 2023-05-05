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
    type=c(rep('Observed', length(Y.test)), rep('Predicted', length(Y.pred)))
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
#one_time_data_clean()

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

n.train <- 2 #floor(0.5*n.total)
b.X.train <- list()
X.train.long <- X.train[[1]]
Y.train.long <- Y.train[[1]]
for (i in 1:n.train) {
  cat(i,' ')
  b.X.train[[i]] <- basis_maker(X.all[[i]])
  
  if (1 != i) {
    X.train.long <- rbind(X.train.long, X.all[[i]])
    Y.train.long <- c(Y.train.long, Y.all[[i]])
  }
}

kernlearn.out <- fit_kernlearn_Xunique(b.X.train, basis_maker, Y.all[1:n.train])

E.beta <- kernlearn.out$E.est
V.beta <- kernlearn.out$V.est

test.itr <- 4
X.test <- X.all[[test.itr]]
post.out <- posterior_test_meanvar_brev(E.beta, V.beta, 
                                        X.train.long, Y.train.long, 
                                        X.test, basis_maker)
p <- plot3D_piped(X.train[[4]], Y.train[[4]], post.out$post.test.mean)
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

