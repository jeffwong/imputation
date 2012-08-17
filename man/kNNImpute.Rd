\name{kNNImpute}
\alias{kNNImpute}
\title{
k-NN Imputation
}
\description{
Imputation using k nearest neighbors
}
\usage{
kNNImpute(x, k, verbose=T)
}
\arguments{
  \item{x}{
Numeric data matrix containing missing values
}
  \item{k}{
The number of neighbors to use for k-NN
}
  \item{verbose}{
Logical.  If true, print updates about the status of kNNImpute
}
}
\details{
Identifies observations, x, that have missing data.  If column j of observation 
x is missing, finds the k nearest observations to x that have non-missing 
values in j.  Among this set, use the mean of j to fill the missing value in x.
Reports warnings if an entire row is missing.
}
\value{
\item{x }{The imputed data matrix}
\item{missing.matrix }{A matrix of logicals displaying the location of the original missing values}
}
\references{
}
\author{
Jeffrey Wong
}
\note{
}

\seealso{
SVDImpute, SVTImpute
}
\examples{
x = matrix(rnorm(30),10,3)
x[2,3] = NA
x[5,1] = NA
x[2,2] = NA
x.imputed = kNNImpute(x,2)
}
\keyword{ k-NN }
\keyword{ impute }
