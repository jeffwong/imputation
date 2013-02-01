\name{cv.lmImpute}
\alias{cv.lmImpute}
\title{CV for lmImpute}
\usage{
  cv.lmImpute(x, ...)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{...}{additional parameters passed to locfit}
}
\description{
  Cross Validation for Locally Weighted Linear Imputation
  Artificially erase some data and run lmImpute to compute
  the RMSE on the subset of x for which data was
  artificially erased.
}

