\name{cv.SVDImpute}
\alias{cv.SVDImpute}
\title{
Cross Validation for SVDImpute
}
\description{
Use cross-validation to determine optimal parameters for kNN, SVD, and SVD imputation
}
\usage{
cv.SVDImpute(x, k.max = floor(ncol(x)/2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data matrix containing missing values.
}
  \item{k.max}{
The maximum rank-k approximation that should be used
}
}
\details{
Erases data at random and attempts SVD imputation to restore the matrix.
The mean relative absolute error is used to determine what value of k gives the
best results.
}
\value{
\item{k }{The optimal k to use determined by cross validation}
\item{mae }{The mean relative absolute error when using a rank-k approximation}
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
cv.kNNImpute, cv.SVTImpute
}
\examples{
x = matrix(rnorm(100),10,10)
x[1,1] = NA
x[5,8] = NA
x[2,3] = NA
x[7,7] = NA
k = cv.SVDImpute(x)$k
SVDImpute(x, k)
}
