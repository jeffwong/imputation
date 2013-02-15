\name{SVDImpute}
\alias{SVDImpute}
\title{SVD Imputation}
\usage{
  SVDImpute(x, k, num.iters = 10, verbose = T)
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k}{the rank-k approximation to use for x}

  \item{num.iters}{the number of times to compute the
  rank-k approximation and impute the missing data}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using the SVD First fill missing values using
  the mean of the column Then, compute a low, rank-k
  approximation of x.  Fill the missing values again from
  the rank-k approximation.  Recompute the rank-k
  approximation with the imputed values and fill again,
  repeating num.iters times
}
\examples{
x = matrix(rnorm(100),10,10)
  x.missing = x > 1
  x[x.missing] = NA
  SVDImpute(x, 3)
}

