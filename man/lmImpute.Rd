\name{lmImpute}
\alias{lmImpute}
\title{Locally Weighted Linear Imputation}
\usage{
  lmImpute(x)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}
}
\description{
  Fill missing values in a column by running a locally
  weighted least squares regression against the row number.
  Good for large data (large number of records)
}

