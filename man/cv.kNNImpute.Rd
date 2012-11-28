\name{cv.kNNImpute}
\alias{cv.kNNImpute}
\title{CV for kNNImpute}
\usage{
  cv.kNNImpute(x, k.max = 5)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k.max}{the largest amount of neighbors to try kNN
  Impute}
}
\description{
  Cross Validation for kNNImpute Artificially erase some
  data and run kNNImpute multiple times, varying k from 1
  to k.max.  For each k, compute the RMSE on the subset of
  x for which data was artificially erased.
}

