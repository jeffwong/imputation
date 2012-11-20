kNNImpute = function(x, k, x.dist = NULL, impute.fn = mean, verbose=T) {
  if(k >= nrow(x))
    stop("k must be less than the number of rows in x")

  prelim = impute.prelim(x)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$x.missing

  if (verbose) print("Computing distance matrix...")
  if (is.null(x.dist) & nrow(x) < 2000) x.dist = as.matrix(dist(x, upper=T))
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
      if (!is.null(x.dist))
        knn.ranks = order(x.dist[rowIndex,neighbor.indices])
      else
        knn.ranks = order(as.vector(pdist(x, rowIndex, neighbor.indices)))
      #identify the row number in the original data matrix of the knn
      knn = neighbor.indices[(knn.ranks[1:k])]
      impute.fn(x[knn,j])
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

cv.kNNImpute = function(x, k.max=5) {
  if(k.max >= nrow(x)) stop("k.max must be less than nrow(x)")

  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.dist = if (nrow(x) <= 2000) dist(x, upper = T) else NULL
  rmse = sapply(1:k.max, function(i) {
    x.imputed = kNNImpute(x.train, i, x.dist, verbose=F)$x
    error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
    sqrt(mean(error^2))
  })
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse)
}
