rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('./package/gp_kernlearn/code/model_choice_ic.R')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get AIC/BIC for Our Models ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
models <- list.files('./package/gp_kernlearn/code/output/gof_bs/models/')
dataset.reps <- 100

output <- data.frame(
  model=rep(models, each=dataset.reps),
  data_gen=rep(NA, dataset.reps*length(models)),
  dataset_num=rep(NA, dataset.reps*length(models)),
  leg_deg=rep(NA, dataset.reps*length(models)),
  llik=rep(NA, dataset.reps*length(models)),
  aic=rep(NA, dataset.reps*length(models)),
  bic=rep(NA, dataset.reps*length(models))
)

for (i in 1:nrow(output)) {
  cat(i, ' ')
  name <- strsplit(output$model[i], '.RData')[[1]]
  output$data_gen[i] <- strsplit(name, '_')[[1]][2]
  
  leg.deg.chr <- strsplit(name, '_')[[1]][length(strsplit(name, '_')[[1]])]
  leg.deg <- as.numeric(substr(leg.deg.chr, 5, nchar(leg.deg.chr)))
  output$leg_deg[i] <- leg.deg
  
  model <- readRDS(paste0('./package/gp_kernlearn/code/output/gof_bs/models/',output$model[i]))

  dataset.itr <- ifelse(0==i%%dataset.reps, dataset.reps, i%%dataset.reps)
  output$dataset_num[i] <- dataset.itr
  n <- length(model$Y.train[[dataset.itr]])
  llik <- mvn_loglikelihood(model$Y.train[[dataset.itr]][,1], mu = rep(0, n), model$K.est)
  output$llik[i] <- llik
  output$aic[i] <- aic(llik, leg.deg)
  output$bic[i] <- bic(llik, leg.deg, n)
  
}

write.csv(output, paste0('./package/gp_kernlearn/code/output/modchoice_ic/', 'output.csv'), row.names = FALSE)
output <- read.csv(paste0('./package/gp_kernlearn/code/output/modchoice_ic/', 'output.csv'))
#output$aic <- output$aic + abs(min(output$aic)) + 1
#output$bic <- output$bic + abs(min(output$bic)) + 1

output <- output[output$data_gen != 'Ylin',]
output$facet_name <- factor(output$data_gen, 
                            levels=c('Ylin1-1', 'Ycubic', 'Yquad1-1-1', 
                                     'Yrbf5-10', 'Yrbf5-05-fz01','Yrbf5-005-fz01'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get AIC/BIC for Our Models (stationary) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source('package/gp_kernlearn/code/stationary.R')

models <- list.files('./package/gp_kernlearn/code/output/gof_bs/models/')
dataset.reps <- 100

output.stationary <- data.frame(
  model=rep(models, each=dataset.reps),
  data_gen=rep(NA, dataset.reps*length(models)),
  dataset_num=rep(NA, dataset.reps*length(models)),
  leg_deg=rep(NA, dataset.reps*length(models)),
  llik_stationary=rep(NA, dataset.reps*length(models)),
  aic_stationary=rep(NA, dataset.reps*length(models)),
  bic_stationary=rep(NA, dataset.reps*length(models))
)

for (i in 1:nrow(output.stationary)) {
  cat(i, ' ')
  name <- strsplit(output.stationary$model[i], '.RData')[[1]]
  output.stationary$data_gen[i] <- strsplit(name, '_')[[1]][2]
  
  leg.deg.chr <- strsplit(name, '_')[[1]][length(strsplit(name, '_')[[1]])]
  leg.deg <- as.numeric(substr(leg.deg.chr, 5, nchar(leg.deg.chr)))
  output.stationary$leg_deg[i] <- leg.deg
  
  model <- readRDS(paste0('./package/gp_kernlearn/code/output/gof_bs/models/',output.stationary$model[i]))
  
  dataset.itr <- ifelse(0==i%%dataset.reps, dataset.reps, i%%dataset.reps)
  output.stationary$dataset_num[i] <- dataset.itr
  n <- length(model$Y.train[[dataset.itr]])
  
  X <- model$b.X.train[,2]
  s.out.mon3 <- make_stationary(model$K.est, X, dist_euc, modfit_monopoly_maker(3))
  llik <- mvn_loglikelihood(model$Y.train[[dataset.itr]][,1], mu = rep(0, n), s.out.mon3$K.stationary)
  output.stationary$llik_stationary[i] <- llik
  output.stationary$aic_stationary[i] <- aic(llik, leg.deg)
  output.stationary$bic_stationary[i] <- bic(llik, leg.deg, n)
}

write.csv(output.stationary, paste0('./package/gp_kernlearn/code/output/modchoice_ic/', 'output_stationary_mon3.csv'), row.names = FALSE)

output.stationary <- read.csv(paste0('./package/gp_kernlearn/code/output/modchoice_ic/', 'output_stationary_mon3.csv'))

output.stationary <- output.stationary[output.stationary$data_gen != 'Ylin',]
output.stationary$facet_name <- factor(output.stationary$data_gen, 
                            levels=c('Ylin1-1', 'Ycubic', 'Yquad1-1-1', 
                                     'Yrbf5-10', 'Yrbf5-05-fz01','Yrbf5-005-fz01'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Fit Basic GP ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gplite)

models <- list.files('./package/gp_kernlearn/code/output/gof_bs/models/')
dataset.reps <- 100

output.gp <- data.frame(
  model=models,
  data_gen=rep(NA, length(models)),
  llik=rep(NA, length(models)),
  aic=rep(NA, length(models)),
  bic=rep(NA, length(models))
)

for (i in 1:nrow(output.gp)) {
  cat(i, ' ')
  rbf.df <- 2
  name <- strsplit(output.gp$model[i], '.RData')[[1]]
  datagen.name <- strsplit(name, '_')[[1]][2]
  output.gp$data_gen[i] <- datagen.name 
  
  if (1 == i ) {
    covered.names <- c(datagen.name)
  } else if (datagen.name %in% covered.names) {
    next
  } else {
    covered.names <- c(covered.names, datagen.name)
  }
  
  model <- readRDS(paste0('./package/gp_kernlearn/code/output/gof_bs/models/',output.gp$model[i]))
  
  
  for (dataset.itr in 1:dataset.reps) {
    gp <- gp_init(cfs = cf_sexp(), lik = lik_gaussian())
    gp <- gp_optim(gp, model$b.X.train[,2], model$Y.train[[dataset.itr]], verbose=FALSE)
    
    if (1 == dataset.itr) {
      K.gp <- gp$fit$C_chol%*%t(gp$fit$C_chol)  
    } else {
      K.gp <- K.gp + gp$fit$C_chol%*%t(gp$fit$C_chol)  
    }
  }
  
  K.gp.avg <- K.gp/dataset.reps
  # TODO image(K.gp.avg) RBF should not work for cubic data, for example
  n <- length(model$Y.train[[dataset.itr]])
  llik <- mvn_loglikelihood(model$Y.train[[dataset.itr]][,1], mu = rep(0, n),K.gp.avg)
  output.gp$llik[i] <- llik
  output.gp$aic[i] <- aic(llik, rbf.df)
  output.gp$bic[i] <- bic(llik, rbf.df, n)
  break
}

output.gp <- output.gp[!is.na(output.gp$llik),]
output.gp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plotting ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
ggplot(data=output.stationary) + 
  geom_line(aes(x=leg_deg, y=aic_stationary, group=dataset_num, color=data_gen), alpha=0.5) +
  facet_wrap(~facet_name, scale='free') +
  theme(legend.position = 'bottom') #+
#  scale_y_continuous(trans='log10')


ggplot(data=output.stationary) + 
  geom_line(aes(x=leg_deg, y=bic_stationary, group=dataset_num, color=data_gen), alpha=0.5) +
  facet_wrap(~facet_name, scale='free') +
  theme(legend.position = 'bottom')
