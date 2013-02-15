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
\examples{
dates = timeSequence(from = '2012-01-01', to = '2012-12-31', by = 'day')
  dimensions = sample(c("A", "B"), 366, replace = TRUE)
  numA = length(which(dimensions == "A")); numB = length(which(dimensions == "B"))
  metrics = matrix(0, 366, 2)
  metrics[which(dimensions == "A"),1] = rnorm(numA, mean=1)
  metrics[which(dimensions == "A"),2] = rnorm(numA, mean=5)
  metrics[which(dimensions == "B"),1] = rnorm(numB, mean=-10)
  metrics[which(dimensions == "B"),2] = rnorm(numB, mean=-5)
  tp = projectDate(as.Date(dates))
  monday.indices = which(tp$weekday == "Monday")
  metrics[sample(monday.indices, 20),] = NA
  cv.tsImpute(as.Date(dates), dimensions, metrics)
}

