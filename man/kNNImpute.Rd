\name{kNNImpute}
\alias{kNNImpute}
\title{kNN Impute}
\usage{
  kNNImpute(x, k, x.dist = NULL, impute.fn = mean, ...,
    verbose = T)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k}{the number of neighbors to use for imputation}

  \item{x.dist}{an optional, pre-computed distance matrix
  to be used for kNN}

  \item{impute.fn}{the imputation function to run on the
  length k vector of values for a missing feature}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using k-nearest neighbors For each record,
  identify missinng features.  For each missing feature
  find the k nearest neighbors which have that feature.
  Impute the missing value using the imputation function on
  the k-length vector of values found from the neighbors
}

