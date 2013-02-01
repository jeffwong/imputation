\name{robustSVDImpute}
\alias{robustSVDImpute}
\title{Robust SVD Imputation}
\usage{
  robustSVDImpute(x, k, alpha = 1/2, max.iters = 10,
    verbose = T)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k}{the rank-k approximation to use for x}

  \item{alpha}{the alpha passed into least trimmed squares}

  \item{max.iters}{the number of times to alternate least
  trimmed squares}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using a robust implementation of the SVD
  Formulate the SVD as a biconvex minimization problem, and
  solve it by alternating least squares.  When computing
  least squares, use the least trimmed squares
  implementation, a robust implementation that trims
  extreme values out
}

