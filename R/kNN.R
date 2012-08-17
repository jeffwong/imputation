kNNImpute = function(x, k, x.dist = NULL, impute.fn = mean, verbose=T) {
  if(k >= nrow(x))
    stop("k must be less than the number of rows in x")
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix) 
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
      nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return (x)
  }
  
  if(verbose)
    print("Computing distance matrix...")
  if (is.null(x.dist)) {
      x.dist = as.matrix(dist(x, upper=T))
  }
  if(verbose)
    print("Distance matrix complete")
  
  missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
    any(i)
  }))
  x.missing = (cbind(1:nrow(x),x))[missing.rows.indices,]
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
      knn.ranks = order(x.dist[rowIndex,neighbor.indices])
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
  if(k.max >= nrow(x))
    stop("k.max must be less than nrow(x)")
  n = nrow(x) * ncol(x)
  missing.matrix = is.na(x)
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, 1/3*length(valid.data))
  x.train = x
  x.train[remove.indices] = NA

  absolute.error = sapply(1:k.max, function(i) {
    x.imputed = kNNImpute(x.train, i, verbose=F)$x
    mae = abs((x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices] )
  })
  mae = apply(absolute.error, 2, function(j) {
    mean(j)
  })
  list(k = which.min(mae), mae = mae[which.min(mae)],
    k.full = 1:k.max, mae.full = mae)
}
