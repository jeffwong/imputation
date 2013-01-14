\name{cv.SVTImpute}
\alias{cv.SVTImpute}
\title{CV for SVTImpute}
\usage{
  cv.SVTImpute(x,
    lambda.range = seq(0, 1, length.out = 101),
    parallel = F)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{lambda.range}{a vector of penalty terms to use in
  the CV}
}
\description{
  Cross Validation for SVT Imputation Artificially erase
  some data and run SVTImpute multiple times, varying
  lambda using lambda.range.  For each lambda, compute the
  RMSE on the subset of x for which data was artificially
  erased.
}

