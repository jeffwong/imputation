#' Locally Weighted Linear Imputation
#'
#' Fill missing values in a column by running a locally weighted least squares
#' regression against the row number.
#' Good for large data (large number of records)
#' @param x a data frame or matrix where each row represents a different record
#' @param ... additional parameters passed to locfit
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   lmImpute(x)
#' @export
lmImpute = function(x, ...) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing

  indices = 1:nrow(x)
  x.imputed = apply(x.missing, 2, function(j) {
    predict(locfit(j[-1] ~ indices, ...), indices)
  })
  x[,prelim$missing.cols.indices] = x.imputed
  
  return (list(x = x,
               missing.matrix = missing.matrix))

}

#' CV for lmImpute
#'
#' Cross Validation for Locally Weighted Linear Imputation
#' Artificially erase some data and run lmImpute to compute the RMSE on 
#' the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param ... additional parameters passed to locfit
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   cv.lmImpute(x)
#' @export
cv.lmImpute = function(x, ...) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.imputed = lmImpute(x.train)$x
  error = (x[remove.indices] - x.imputed[remove.indices])
  nerror = error / x[remove.indices]
  rmse = sqrt(mean(error^2))
  nrmse = sqrt(mean(nerror^2))
  
  list(imputation = x.imputed, rmse = rmse, nrmse = nrmse)
}
