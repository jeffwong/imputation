\name{lmImpute}
\alias{lmImpute}
\title{Locally Weighted Linear Imputation}
\usage{
  lmImpute(x, ...)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{...}{additional parameters passed to locfit}
}
\description{
  Fill missing values in a column by running a locally
  weighted least squares regression against the row number.
  Good for large data (large number of records)
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  lmImpute(x)
}

