\name{matern_kernel_euc}
\alias{matern_kernel_euc}
\title{Create Matern Kernel Matrix}
\description{
  Evalute the Matern covariance function at the specified values.
}
\usage{
  matern_kernel_euc(X, sig.0, rho, nu, fuzz)
}
\arguments{
  \item{X}{number of independent samples/datasets to generate}
  \item{sig.0}{overall variance modifier}
  \item{rho}{}
  \item{nu}{}
  \item{fuzz}{variance of in}
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
}
\keyword{ultility, util}
