#' kNN Impute
#'
#' Imputation using k-nearest neighbors
#' For each record, identify missinng features.  For each missing feature
#' find the k nearest neighbors which have that feature.  Impute the missing
#' value using the imputation function on the k-length vector of values
#' found from the neighbors
#' @param x a data frame or matrix where each row represents a different record
#' @param k the number of neighbors to use for imputation
#' @param x.dist an optional, pre-computed distance matrix to be used for kNN
#' @param impute.fn the imputation function to run on the length k vector of values for
#'   a missing feature
#' @param verbose if TRUE print status updates
#' @export
kNNImpute = function(x, k, x.dist = NULL, impute.fn = mean, ..., verbose=T) {
  if(k >= nrow(x))
    stop("k must be less than the number of rows in x")

  prelim = impute.prelim(x)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing
  missing.rows.indices = prelim$missing.rows.indices

  if (verbose) print("Computing distance matrix...")
  if (is.null(x.dist)) x.dist = dist(x)
  if (verbose) print("Distance matrix complete")
  
  x.missing.imputed = t(apply(x.missing, 1, function(i) {
    rowIndex = i[1]
    i.original = i[-1]
    if(verbose) print(paste("Imputing row", rowIndex,sep=" "))
    missing.cols = which(missing.matrix[rowIndex,])
    if(length(missing.cols) == ncol(x))
      warning( paste("Row",rowIndex,"is completely missing",sep=" ") )
    imputed.values = sapply(missing.cols, function(j) {
      #find neighbors that have data on the jth column
      neighbor.indices = which(!missing.matrix[,j])
      #lookup the distance to these neighbors
      #order the neighbors to find the closest ones
      if (!is.null(x.dist)) {
        indices.1d = .dist.2dto1d(rowIndex, neighbor.indices, nrow(x))
        knn.ranks = order(x.dist[indices.1d])
      }
      else
        knn.ranks = order(pdist(x, indices.A = rowIndex,
                          indices.B = neighbor.indices)@dist)
      #identify the row number in the original data matrix of the knn
      knn = neighbor.indices[(knn.ranks[1:k])]
      impute.fn(x[knn,j], ...)
    })
    i.original[missing.cols] = imputed.values
    i.original
  }))
  x[missing.rows.indices,] = x.missing.imputed

  missing.matrix2 = is.na(x)
  x[missing.matrix2] = 0

  return (list(
    x=x,
    missing.matrix=missing.matrix
  ))
}

#' CV for kNNImpute
#'
#' Cross Validation for kNNImpute
#' Artificially erase some data and run kNNImpute multiple times,
#' varying k from 1 to k.max.  For each k, compute the RMSE on the subset of x
#' for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param k.max the largest amount of neighbors to try kNN Impute
#' @export
cv.kNNImpute = function(x, k.max=5) {
  if(k.max >= nrow(x)) stop("k.max must be less than nrow(x)")

  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.dist = dist(x)
  rmse = sapply(1:k.max, function(i) {
    x.imputed = kNNImpute(x.train, i, x.dist, verbose=F)$x
    error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
    sqrt(mean(error^2))
  })
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse)
}

#' 2D indices to 1D indices
#'
#' Helper function to convert 2D indices to 1D indices.
#' The return value of the function dist does not by default return a
#' 2D object, instead it returns an array.  When wanting to access an
#' element at the i,jth position of the distance matrix, this function
#' converts the 2D index to a 1D index that can be used on the distance array
.dist.2dto1d = function(i,j,n) {
  ret = rep(0, length(j))
  j.larger.indices = which(j > i)
  j.smaller.indices = which(j < i)
  if (length(j.larger.indices) > 0) {
    j.larger = j[j.larger.indices]
    ret[1:length(j.larger)] = (j.larger-1)*n - j.larger^2/2 + i - 1
  }
  if (length(j.smaller.indices) > 0) {
    j.smaller = j[j.smaller.indices]
    ret[(length(j.larger.indices)+1) : length(j)] = (i-1)*n - i^2/2 + j.smaller - 1
  }
  ret
}
