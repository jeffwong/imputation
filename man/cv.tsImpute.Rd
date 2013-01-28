\name{cv.tsImpute}
\alias{cv.tsImpute}
\title{CV for tsImpute}
\usage{
  cv.tsImpute(time, dimension, metric, ...)
}
\arguments{
  \item{time}{a vector of dates or datetime objects}

  \item{dimension}{a data frame of exogenous predictor
  variables}

  \item{metric}{a matrix where each column represents a
  time series}

  \item{...}{extra parameters to be passed to tsImpute}
}
\description{
  Cross Validation for Time Series Imputation Artificially
  erase some data and run gbmImpute.  Compute the RMSE on
  the subset of x for which data was artificially erased.
}

