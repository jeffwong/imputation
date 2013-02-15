\name{cv.SVTApproxImpute}
\alias{cv.SVTApproxImpute}
\title{CV for SVTApproxImpute}
\usage{
  cv.SVTApproxImpute(x,
    lambda.range = seq(0, 1, length.out = 101),
    parallel = F)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{lambda.range}{a vector of penalty terms to use in
  the CV}

  \item{parallel}{runs each run for lambda in lambda.range
  in parallel.  Requires a parallel backend to be
  registered}
}
\description{
  Cross Validation for SVT Imputation Artificially erase
  some data and run SVTImpute multiple times, varying
  lambda using lambda.range.  For each lambda, compute the
  RMSE on the subset of x for which data was artificially
  erased.
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  cv.SVTApproxImpute(x)
}

