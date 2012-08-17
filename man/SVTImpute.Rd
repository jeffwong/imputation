\name{SVTImpute}
\alias{SVTImpute}
\title{
SVT Imputation
}
\description{
Imputation using Singular Value Thresholding
}
\usage{
SVTImpute(x, lambda, verbose = F)
}
\arguments{
  \item{x}{
Numeric matrix containing some missing data.  Each row is an observation
and the columns are features of the data
}
  \item{lambda}{
A tuning parameter used for thresholding.  See details
}
  \item{verbose}{
Logical.  If true, print updates about the status of SVTImpute
}
}
\details{
Uses Singular Value Thresholding for Imputation.  Computes the SVD of x and
an augmented singular value matrix by subtracting lambda from each singular 
value of x.  If this results in any values less than 0, set these values equal
to 0.  Reconstruct x by multiplying through the SVD using the augmented singular
value matrix.
}
\value{
\item{x }{The imputed data matrix}
\item{missing.matrix }{A matrix of logicals displaying the location of the original missing values}
}
\references{
Emmanuel Candes, Stanford University
}
\author{
Jeffrey Wong
}
\note{

}

\seealso{
kNNImpute, SVDImpute
}
\examples{
x = matrix(rnorm(30),10,3)
x[2,3] = NA
x[5,1] = NA
x[2,2] = NA
x.imputed = SVTImpute(x,0.1)
}
