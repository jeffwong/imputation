#' Mean Imputation
#'
#' Fill missing values in a column with the mean of the column
#' @param x a data frame or matrix where each row represents a different record
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   meanImpute(x)
#' @export
meanImpute = function(x) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing

  x.imputed = apply(x.missing, 2, function(j) {
      bad.indices = which(is.na(j))
      j[bad.indices] = mean(j[-1], na.rm = T)
      j[-1]
  })
  x[,prelim$missing.cols.indices] = x.imputed
  
  return (list(x = x,
               missing.matrix = missing.matrix))

}

#' CV for meanImpute
#'
#' Cross Validation for mean Imputation
#' Artificially erase some data and run meanImpute to compute the RMSE on 
#' the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   cv.meanImpute(x)
#' @export
cv.meanImpute = function(x) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.imputed = meanImpute(x.train)$x
  error = (x[remove.indices] - x.imputed[remove.indices])
  nerror = error / x[remove.indices]
  rmse = sqrt(mean(error^2))
  nrmse = sqrt(mean(nerror^2))
  
  list(imputation = x.imputed, rmse = rmse, nrmse = nrmse)
}
