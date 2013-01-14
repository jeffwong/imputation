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
}
\description{
  Cross Validation for SVD Imputation Artificially erase
  some data and run SVDImpute multiple times, varying k
  from 1 to k.max.  For each k, compute the RMSE on the
  subset of x for which data was artificially erased.
}

