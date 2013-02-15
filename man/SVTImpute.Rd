\name{SVTImpute}
\alias{SVTImpute}
\title{SVT Imputation}
\usage{
  SVTImpute(x, lambda, stepsize, threshold = 0.001,
    max.iters = 10, verbose = F)
}
\arguments{
  \item{x}{a data frame or matrix with size n1 x n2 where
  each row represents a different record}

  \item{lambda}{the penalty on the singular values}

  \item{stepsize}{optional.  If not provided, uses 1.2 *
  (n1 * n2) / (number of missing elements)}

  \item{threshold}{convergence threshold}

  \item{max.iters}{maximum number of iterations.  Note that
  each iteration will require computing an SVD}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using Singular Value Thresholding a la Cai,
  Candes, Shen.
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  SVTImpute(x, 3)
}
\references{
  A Singular Value Thresholding Algorithm for Matrix
  Completion. Cai, Candes, Shen.
}

