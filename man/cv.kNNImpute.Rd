\name{cv.kNNImpute}
\alias{cv.kNNImpute}
\title{
Cross Validation for kNNImpute
}
\description{
Use cross-validation to determine optimal parameters for kNN, SVD, and SVD imputation
}
\usage{
cv.kNNImpute(x, k.max = 5)
}
\arguments{
  \item{x}{
Data matrix containing missing values.
}
  \item{k.max}{
The maximum number of neighbors that should be attempted
}
}
\details{
Erases data at random and attempts kNN imputation to restore the matrix.
The mean relative absolute error is used to determine what value of k gives the
best results.
}
\value{
\item{k }{The optimal number of neighbors to use determined by cross validation}
\item{mae }{The mean relative absolute error when using k neighbors}
\item{k.full }{All values of k that were attempted}
\item{mae.full }{All mean relative absolute errors that were inspected}
}
\references{

}
\author{
Jeffrey Wong
}
\note{

}

\seealso{
cv.SVDImpute, cv.SVDImpute
}
\examples{
x = matrix(rnorm(100),10,10)
x[1,1] = NA
x[5,8] = NA
x[2,3] = NA
x[7,7] = NA
k = cv.kNNImpute(x)$k
kNNImpute(x, k)
}

