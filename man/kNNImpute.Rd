\name{kNNImpute}
\alias{kNNImpute}
\title{kNN Impute}
\usage{
  kNNImpute(x, k, x.dist = NULL, impute.fn, ...,
    verbose = T)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k}{the number of neighbors to use for imputation}

  \item{x.dist}{an optional, pre-computed distance matrix
  to be used for kNN}

  \item{impute.fn}{the imputation function to run on the
  length k vector of values for a missing feature.
  Defaults to a weighted mean of the neighboring values
  weighted by the distance of the neighbors}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using k-nearest neighbors. For each record,
  identify missinng features.  For each missing feature
  find the k nearest neighbors which have that feature.
  Impute the missing value using the imputation function on
  the k-length vector of values found from the neighbors
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 2
  x[x.missing] = NA
  kNNImpute(x, 3)
}

