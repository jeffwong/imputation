\name{meanImpute}
\alias{meanImpute}
\title{Mean Imputation}
\usage{
  meanImpute(x)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}
}
\description{
  Fill missing values in a column with the mean of the
  column
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  meanImpute(x)
}

