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
  The default impute.fn weighs the k values by their
  respective distances.  First the smallest k distances are
  extracted into the variable smallest.distances Then, the
  corresponding values are extracted to knn.values.
  Finally, knn.weights normalizes the distances by the max
  distance, and are subtracted by 1.  The result is the
  weighted mean of the values of the nearest neighbors and
  their weight based on their distance.  It is implemented
  as follows: \code{impute.fn = function(values, distances,
  k) { ranks = order(distances) smallest.distances =
  distances[ranks][1:k] #values corresponding to smallest
  distances knn.values = values[ranks][1:k] knn.weights = 1
  - (smallest.distances / max(distances))
  weighted.mean(knn.values, knn.weights) }}
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  kNNImpute(x, 3)
}
\references{
  Missing value estimation methods for DNA microarrays.
  Troyanskaya et al.
}

