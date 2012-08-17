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
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix)
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
      nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return (x)
  }
  missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
    any(i)
  }))
  x.missing = (rbind(1:ncol(x), x))[,missing.cols.indices]
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
  n = nrow(x) * ncol(x)
  missing.matrix = is.na(x)
  if(sum(missing.matrix) == 0) stop("No missing data")
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, 1/3*length(valid.data))
  x.train = x
  x.train[remove.indices] = NA

  absolute.error = sapply(1:k.max, function(i) {
    x.imputed = SVDImpute(x.train, i, verbose=F)$x
    mae = abs( (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices] )
  })
  mae = apply(absolute.error, 2, function(j) {
    mean(j)
  })
  list(k = which.min(mae), mae = mae[which.min(mae)],
    k.full = 1:k.max, mae.full = mae)
}

SVTImpute = function(x, lambda, verbose=F) {
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix)
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
      nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return (x)
  }
  missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
    any(i)
  }))
  x.missing = (rbind(1:ncol(x), x))[,missing.cols.indices]
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

cv.SVTImpute = function(x, lambda.range = seq(0,1,length.out=101)) {
  n = nrow(x) * ncol(x)
  missing.matrix = is.na(x)
  if(sum(missing.matrix) == 0) stop("No missing data")
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, 1/3*length(valid.data))
  x.train = x
  x.train[remove.indices] = NA

  absolute.error = sapply(lambda.range, function(i) {
    x.imputed = SVTImpute(x.train, i, verbose=F)$x
    mae = abs( (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices] )
  })
  mae = apply(absolute.error, 2, function(j) {
    mean(j)
  })
  list(lambda = lambda.range[which.min(mae)], mae = mae[which.min(mae)],
    lambda.full = lambda.range, mae.full = mae)
}

