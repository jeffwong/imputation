.rankKapprox = function(x, k) {
  x.svd = svd(x, nu=k, nv=k)
  x.svd$u %*% diag(x.svd$d[1:k],nrow=k,ncol=k) %*% t(x.svd$v)
}

#' SVD Imputation
#'
#' Imputation using the SVD
#' First fill missing values using the mean of the column
#' Then, compute a low, rank-k approximation of x.  Fill the
#' missing values again from the rank-k approximation.  Recompute
#' the rank-k approximation with the imputed values and fill again,
#' repeating num.iters times
#' @param x a data frame or matrix where each row represents a different record
#' @param k the rank-k approximation to use for x
#' @param num.iters the number of times to compute the rank-k approximation
#'   and impute the missing data
#' @param verbose if TRUE print status updates
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   SVDImpute(x, 3)
#' @export
SVDImpute = function(x, k, num.iters = 10, verbose=T) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing
  missing.cols.indices = prelim$missing.cols.indices

  #First initialize missing values with mean
  x.missing.imputed = apply(x.missing, 2, function(j) {
    colIndex = j[1]
    j.original = j[-1]
    missing.rows = which(missing.matrix[,colIndex])
    if(length(missing.rows) == nrow(x))
      warning( paste("Column",colIndex,"is completely missing",sep=" ") )
    j.original[missing.rows] = mean(j.original[-missing.rows])
    j.original
  })
  #replace columns with missing values with x.missing.imputed
  x[,missing.cols.indices] = x.missing.imputed
  #Fill anything that is still NA with 0
  missing.matrix2 = is.na(x)
  x[missing.matrix2] = 0
  for(i in 1:num.iters) {
    if(verbose) print(paste("Running iteration", i, sep=" "))
    x.svd = .rankKapprox(x, k)
    x[missing.matrix] = x.svd[missing.matrix]
  }
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}

#' CV for SVDImpute
#'
#' Cross Validation for SVD Imputation
#' Artificially erase some data and run SVDImpute multiple times,
#' varying k from 1 to k.max.  For each k, compute the RMSE on the subset of x
#' for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param k.max the largest rank used to approximate x
#' @param parallel runs each run for k = 1 to k = k.max in parallel.  Requires
#'   a parallel backend to be registered
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   cv.SVDImpute(x)
#' @export
cv.SVDImpute = function(x, k.max=floor(ncol(x)/2), parallel = F) {
  if(k.max > ncol(x)) {
    stop("Rank-k approximation cannot exceed the number
      of columns of x")
  }

  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  if (parallel) {
    if (!require(foreach)) stop("R package foreach is required for parallel execution, as well
                                 as a registered parallel backend")
    rmse = foreach (i = 1:k.max, .combine = c, .packages = c('imputation')) %dopar% {
      x.imputed = SVDImpute(x.train, i, verbose=F)$x
      error = (x.imputed[remove.indices] - x[remove.indices])
      nerror = error / x[remove.indices]
      list(nrmse = sqrt(mean(nerror^2)), rmse = sqrt(mean(error^2)))
    }
  }
  else {
    rmse = sapply(1:k.max, function(i) {
      x.imputed = SVDImpute(x.train, i, verbose=F)$x
      error = (x.imputed[remove.indices] - x[remove.indices])
      nerror = error / x[remove.indices]
      list(nrmse = sqrt(mean(nerror^2)), rmse = sqrt(mean(error^2))) 
    })
  }
  nrmse = unlist(rmse[1,]); rmse = unlist(rmse[2,])
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse, nrmse.full = nrmse)
}
