\name{cv.persistenceImpute}
\alias{cv.persistenceImpute}
\title{CV for persistenceImpute}
\usage{
  cv.persistenceImpute(x)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different period and columns represent different times
  of that period}
}
\description{
  Cross Validation for persistence Imputation Artificially
  erase some data and run persistenceImpute.  Compute the
  RMSE on the subset of x for which data was artificially
  erased.
}

