\name{cv.gbmImpute}
\alias{cv.gbmImpute}
\title{CV for gbmImpute}
\usage{
  cv.gbmImpute(x, ...)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{...}{extra parameters to be passed to gbmImpute}
}
\description{
  Cross Validation for GBM Imputation Artificially erase
  some data and run gbmImpute.  Compute the RMSE on the
  subset of x for which data was artificially erased.
}
\examples{
x = matrix(rnorm(10000),1000,10)
  x.missing = x > 2
  x[x.missing] = NA
  cv.gbmImpute(x)
}

