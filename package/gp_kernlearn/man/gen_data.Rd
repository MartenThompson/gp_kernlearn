\name{gen_data}
\alias{gen_data}
\title{Generate Gaussian Process Samples (equiv. MV Normal) From Kernel Matrix}
\description{
  The provided kernel is used to generate mean-zero Gaussian process samples. Equivalently, it produces multivariate Normal samples each with mean zero.
}
\usage{
  gen_data(n.train.datasets, Kern)
}
\arguments{
  \item{n.train.datasets}{number of independent samples/datasets to generate}
  \item{Kern}{covariance matrix}
}
\value{
  A list containing the simulated responses.
}
\source{
  Marten Thompson thom7058@umn.edu
}
\examples{
  set.seed(212)
  R <- 20
  n <- 50
  X <- matrix(seq(-1,1, length.out=n), ncol=1)
  sig.0 <- 1
  sig.1 <- 1
  sig.2 <- 1
  fuzz <- 1e-5
  Sigma.true <- matern_kernel_euc(X, sig.0, sig.1, sig.2, fuzz)
  Y.train <- gen_data(R, Sigma.true)
}
\keyword{ultility, util}
