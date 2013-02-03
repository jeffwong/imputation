\name{cv.SVDImpute}
\alias{cv.SVDImpute}
\title{CV for SVDImpute}
\usage{
  cv.SVDImpute(x, k.max = floor(ncol(x)/2), parallel = F)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k.max}{the largest rank used to approximate x}

  \item{parallel}{runs each run for k = 1 to k = k.max in
  parallel.  Requires a parallel backend to be registered}
}
\description{
  Cross Validation for SVD Imputation Artificially erase
  some data and run SVDImpute multiple times, varying k
  from 1 to k.max.  For each k, compute the RMSE on the
  subset of x for which data was artificially erased.
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  cv.SVDImpute(x)
}

