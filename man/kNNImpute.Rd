\name{kNNImpute}
\alias{kNNImpute}
\title{kNN Impute}
\usage{
  kNNImpute(x, k, x.dist = NULL, impute.fn, verbose = T)
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
  the k-length vector of values found from the neighbors.
}
\details{
  The default impute.fn is weighs the k values by their
  respective distances.  It is implemented as follows:
  \code{impute.fn = function(values, distances, k) { ranks
  = order(distances) knn.values = values[order(ranks)][1:k]
  knn.weights = 1 - (distances[ranks] / max(distances))
  [1:k] weighted.mean(knn.values, knn.weights) }}
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 2
  x[x.missing] = NA
  kNNImpute(x, 3)
}
\references{
  Missing value estimation methods for DNA microarrays.
  Troyanskaya et al.
}

