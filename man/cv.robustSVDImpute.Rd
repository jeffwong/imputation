\name{cv.robustSVDImpute}
\alias{cv.robustSVDImpute}
\title{CV for robustSVDImpute}
\usage{
  cv.robustSVDImpute(x, k.max = floor(ncol(x)/2))
}
\arguments{
  \item{x}{a data frame or matrix where each row represents
  a different record}

  \item{k.max}{the largest rank used to approximate x}
}
\description{
  Cross Validation for robust SVD Imputation Artificially
  erase some data and run robustSVDImpute multiple times,
  varying k from 1 to k.max.  For each k, compute the RMSE
  on the subset of x for which data was artificially
  erased.
}

