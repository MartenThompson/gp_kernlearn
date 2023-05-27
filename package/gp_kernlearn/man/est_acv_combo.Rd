\name{est_acv_combo}
\alias{est_acv_combo}
\title{Estimate Autocovariance using Within and Between Sample Variation}
\description{
  This method combines the within-sample ACV with the ACV between the mean responses of all samples. Doing so produces a unique estimate of ACV for each sample \eqn{r=1,\dots,R} and each MCMC sample \eqn{1,\dots, N} for a total of \eqn{R*N} estimates. 
  
  Within-sample ACV is found calling \code{est_acv_within} which performs Bayesian regression.
}
\usage{
  est_acv_combo(Y.list, B.0, n.mcmc, lag, type='covariance')
}
\arguments{
  \item{Y.list}{list of responses across samples, each list element corresponding to a different sample}
  \item{B.0}{design matrix in desired basis}
  \item{n.mcmc}{number of MCMC steps to take when performing initial regression}
  \item{lag}{maximum ACV lag to consider, must be \eqn{\le R, N}}
  \item{type}{ACV method, one of 'correlation' or 'covariance'}
}
\value{
  A list containing the estimated response and estimated ACV values; also, the ACV results from just within-sample analysis.
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
  
  J <- 4
  basis_maker <- make_legendre1D_basis_maker(J-1)
  B.0 <- basis_maker(X)
  
  lag <- 15
  acv.cov.within <- est_acv_combo(Y.train, B.0, 1e2, lag, 'covariance')
}
\keyword{acv}
