\name{mmmImpute}
\alias{mmmImpute}
\title{Mean-Mode-Median Imputation}
\usage{
  mmmImpute(x)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}
}
\description{
  Fill missing values in a column with the mean, mode, or
  median of the column depending on the data type
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  mmmImpute(x)
}

