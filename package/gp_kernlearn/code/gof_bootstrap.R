





bootstrap_stat_dist <- function(statistic, model_gen, n.boot) {
  dummy <- model_gen()
  d <- length(statistic(dummy))
  
  stats <- matrix(NA, n.boot, d)
  
  for (i in 1:n.boot) {
    stats[i,] <- statistic(model_gen())
  }
  
  return(stats)
}

testing <- FALSE
if(testing) {
  m <- function() {
    rnorm(100)
  }
  stat <- function(data) {
    return(c(mean(data), sd(data)))
  }
  
  bs <- bootstrap_stat_dist(stat, m, 100)
  cat(mean(bs$stats[,1]), mean(bs$stats[,2]))
}