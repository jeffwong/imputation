impute.prelim = function(x, byrow = T, verbose=F) {
  missing.matrix = is.na(x)
  numMissing = sum(missing.matrix) 
  if(verbose) {
    print(paste("imputing on", numMissing, "missing values with matrix size",
      nrow(x)*ncol(x), sep=" "))
  }
  if(numMissing == 0) {
    return ( list (missing.matrix = missing.matrix,
                   numMissing = numMissing,
                   missing.rows.indices = NULL,
                   missing.cols.indices = NULL,
                   x.missing = NULL) )
  }

  missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
    any(i)
  }))
  missing.cols.indices = which(apply(missing.matrix, 2, function(i) {
    any(i)
  }))
  if (byrow) x.missing = cbind(1:nrow(x),x)[missing.rows.indices,]
  else x.missing = rbind(1:ncol(x),x)[,missing.cols.indices]

  return ( list (missing.matrix = missing.matrix,
                 numMissing = numMissing,
                 missing.rows.indices = missing.rows.indices,
                 missing.cols.indices = missing.cols.indices,
                 x.missing = x.missing) )
}

cv.impute.prelim = function(x) {
  n = nrow(x) * ncol(x)
  missing.matrix = is.na(x)
  valid.data = which(!missing.matrix)

  remove.indices = sample(valid.data, 1/3*length(valid.data))
  x.train = x; x.train[remove.indices] = NA

  return (list(remove.indices = remove.indices,
               x.train = x.train))
}
