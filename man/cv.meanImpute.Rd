\name{cv.meanImpute}
\alias{cv.meanImpute}
\title{CV for meanImpute}
\usage{
  cv.meanImpute(x)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}
}
\description{
  Cross Validation for mean Imputation Artificially erase
  some data and run meanImpute to compute the RMSE on the
  subset of x for which data was artificially erased.
}

