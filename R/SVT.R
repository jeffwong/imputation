#' SVT Imputation
#'
#' Imputation using Singular Value Thresholding
#' First fill missing values using the mean of the column.
#' Then, compute the SVD of the matrix, and subtract lambda
#' from each of the singular values, thresholding at 0.  Impute
#' by multiplying back out the augmented SVD
#' @param x a data frame or matrix where each row represents a different record
#' @param lambda the penalty on the singular values
#' @param verbose if TRUE print status updates
#' @export
SVTImpute = function(x, lambda, verbose=F) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing
  missing.cols.indices = prelim$missing.cols.indices

  x.missing.imputed = apply(x.missing, 2, function(j) {
    colIndex = j[1]
    j.original = j[-1]
    missing.rows = which(missing.matrix[,colIndex])
    if(length(missing.rows) == nrow(x))
      warning( paste("Column",colIndex,"is completely missing",sep=" ") )
    j.original[missing.rows] = mean(j.original[-missing.rows])
    j.original
  })
  x[,missing.cols.indices] = x.missing.imputed

  #Things that could not be imputed in the initial round should be set to 0
  missing.matrix2 = is.na(x)
  x[missing.matrix2] = 0
   
  x.svd = svd(x)
  lambda.indices = which(x.svd$d < lambda)
  if(length(lambda.indices) > 0) {
    d.augmented = c(x.svd$d[-lambda.indices] - lambda, 
      rep(0, length(lambda.indices)))
  } else {
    d.augmented = x.svd$d - lambda
  }
  x[missing.matrix] = (x.svd$u %*% diag(d.augmented) %*% t(x.svd$v))[missing.matrix]
  return( list(
    x=x,
    missing.matrix = missing.matrix
  ))
}

#' CV for SVTImpute
#'
#' Cross Validation for SVT Imputation
#' Artificially erase some data and run SVTImpute multiple times,
#' varying lambda using lambda.range.  For each lambda, compute the
#' RMSE on the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param lambda.range a vector of penalty terms to use in the CV
#' @param parallel runs each run for lambda in lambda.range in parallel.  Requires
#'   a parallel backend to be registered
#' @export
cv.SVTImpute = function(x, lambda.range = seq(0,1,length.out=101), parallel = F) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  if (parallel) {
    rmse = foreach (i=lambda.range, .combine = c, .packages = c('imputation')) %dopar% {
      x.imputed = SVTImpute(x.train, i, verbose=F)$x
      error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
      sqrt(mean(error^2))
    }
  }
  else {
    rmse = sapply(lambda.range, function(i) {
      x.imputed = SVTImpute(x.train, i, verbose=F)$x
      error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
      sqrt(mean(error^2))
    })
  }
  list(lambda = lambda.range[which.min(rmse)], rmse = rmse[which.min(rmse)],
       lambda.full = lambda.range, rmse.full = rmse)
}
