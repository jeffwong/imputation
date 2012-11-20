#use rank k svd approximation of x
#if x is large, set gpu=T to use CUDA accelerated svd function
.rankKapprox = function(x, k, gpu) {
  if(gpu) {
    x.svd = gpuSvd(x, nu=nrow(x),nv=ncol(x))
  }
  else {
    x.svd = svd(x, nu=k, nv=k)
  }
  x.svd$u %*% diag(x.svd$d[1:k],nrow=k,ncol=k) %*% t(x.svd$v)
}

SVDImpute = function(x, k, num.iters = 10, gpu=F, verbose=T) {
  if(gpu) {
    stop("no gpu support yet")
  }

  prelim = impute.prelim(x)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  x.missing = prelim$missing.matrix

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
    x.svd = .rankKapprox(x, k, gpu)
    x[missing.matrix] = x.svd[missing.matrix]
  }
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}

cv.SVDImpute = function(x, k.max=floor(ncol(x)/2)) {
  if(k.max > ncol(x)) {
    stop("Rank-k approximation cannot exceed the number
      of columns of x")
  }

  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  rmse = sapply(1:k.max, function(i) {
    x.imputed = SVDImpute(x.train, i, verbose=F)$x
    error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
    sqrt(mean(error^2))
  })
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse)
}
