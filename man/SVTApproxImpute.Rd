\name{SVTApproxImpute}
\alias{SVTApproxImpute}
\title{Approximate SVT Imputation}
\usage{
  SVTApproxImpute(x, lambda, verbose = F)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{lambda}{the penalty on the singular values}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using Singular Value Thresholding First fill
  missing values using the mean of the column. Then,
  compute the SVD of the matrix, and subtract lambda from
  each of the singular values, thresholding at 0.  Impute
  by multiplying back out the augmented SVD
}

