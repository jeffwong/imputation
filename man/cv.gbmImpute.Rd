\name{cv.gbmImpute}
\alias{cv.gbmImpute}
\title{CV for gbmImpute}
\usage{
  cv.gbmImpute(date, dimension, metric, ...)

  cv.gbmImpute(date, dimension, metric, ...)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{...}{extra parameters to be passed to gbmImpute}

  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{...}{extra parameters to be passed to gbmImpute}
}
\description{
  Cross Validation for GBM Imputation Artificially erase
  some data and run gbmImpute.  Compute the RMSE on the
  subset of x for which data was artificially erased.

  Cross Validation for Time Series Imputation Artificially
  erase some data and run gbmImpute.  Compute the RMSE on
  the subset of x for which data was artificially erased.
}

