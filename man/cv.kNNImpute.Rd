\name{cv.kNNImpute}
\alias{cv.kNNImpute}
\title{CV for kNNImpute}
\usage{
  cv.kNNImpute(x, k.max = 5, parallel = F)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k.max}{the largest amount of neighbors to try kNN
  Impute}

  \item{parallel}{runs each run for k = 1 to k = k.max in
  parallel.  Requires a parallel backend to be registered}
}
\description{
  Cross Validation for kNNImpute Artificially erase some
  data and run kNNImpute multiple times, varying k from 1
  to k.max.  For each k, compute the RMSE on the subset of
  x for which data was artificially erased.
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  cv.kNNImpute(x)
}

