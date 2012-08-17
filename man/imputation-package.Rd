\name{imputation-package}
\alias{imputation-package}
\alias{imputation}
\docType{package}
\title{
imputation
}
\description{
Provides methods for missing data imputation.  Currently, iterative
SVD imputation, SVT, and kNN are supported.  SVD imputation imputes missing 
values using the mean, then successively calculates a low rank-k approximation 
for the data.  SVT is similar using Singular Value Thresholding.  kNN finds the 
k nearest neighbors and fills missing data by taking the mean of the missing 
variable among the k neighbors.
}
\details{
\tabular{ll}{
Package: \tab imputation\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2011-06-18\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

Very few functions:
\itemize{
  \item SVDImpute - SVD Imputation
  \item cv.SVDImpute - Cross validation to determine optimal k for SVD Imputation
  \item SVTImpute - Singular Value Thresholding Imputation
  \item cv.SVTImpute - Cross validation to determine optimal lambda for SVT Imputation
  \item kNNImpute - k Nearest Neighbors Imputation
  \item cv.kNNImpute - Cross validation to determine optimal k for kNN Imputation
}

}

\author{
Jeffrey Wong
<jeffwong@stanford.edu>
}
\references{
}
\keyword{ package }
\seealso{
}
\examples{
}
