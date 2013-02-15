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
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  cv.meanImpute(x)
}

