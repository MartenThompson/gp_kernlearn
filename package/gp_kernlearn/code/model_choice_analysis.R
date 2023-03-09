rm(list=ls())

setwd('~/Git/gp_kernlearn/')
source('./package/gp_kernlearn/code/model_choice_ic.R')

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

#output$aic <- output$aic + abs(min(output$aic)) + 1
#output$bic <- output$bic + abs(min(output$bic)) + 1

library(ggplot2)
ggplot(data=output) + 
  geom_line(aes(x=leg_deg, y=aic, group=dataset_num, color=data_gen), alpha=0.5) +
  facet_wrap(~data_gen, scale='free')#+
#  scale_y_continuous(trans='log10')


ggplot(data=output) + 
  geom_line(aes(x=leg_deg, y=bic, color=data_gen)) +
  facet_wrap(~data_gen, scale='free')#+
#  scale_y_continuous(trans='log10')
