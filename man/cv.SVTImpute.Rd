\name{cv.SVTImpute}
\alias{cv.SVTImpute}
\title{
Cross Validation for SVTImpute
}
\description{
Use cross-validation to determine optimal parameters for kNN, SVD, and SVD imputation
}
\usage{
cv.SVTImpute(x, lambda.range = seq(0, 1, length.out = 101))
}
\arguments{
  \item{x}{
Data matrix containing missing values.
}
  \item{lambda.range}{
A vector representing the range of lambdas that should be tested
}
}
\details{
Erases data at random and attempts SVT imputation to restore the matrix.
The mean relative absolute error is used to determine what value of lambda gives the
best results.

lambda represents the shrinkage parameter when computing the SVD.  The data matrix
is factored into X = UDV', and the singular values are shrunk by computing max(0, s-lambda)
where s is a singular value of X.  The matrix is restored by multiplying U Dshrink V'
}
\value{
\item{lambda }{The optimal lambda to use determined by cross validation}
\item{mae }{The mean relative absolute error when using lambda}
\item{lambda.full }{All values of lambda that were attempted}
\item{mae.full }{All mean relative absolute errors that were inspected}
}
\references{
http://www-stat.stanford.edu/~candes/papers/SVT.pdf
}
\author{
Jeffrey Wong
}
\note{
}

\seealso{
cv.kNNImpute, cv.SVDImpute
}
\examples{
x = matrix(rnorm(100),10,10)
x[1,1] = NA
x[5,8] = NA
x[2,3] = NA
x[7,7] = NA
lambda = cv.SVTImpute(x)$lambda
SVTImpute(x, lambda)
}
