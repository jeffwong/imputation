impute.prelim = function(x) {
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
  x.missing = (cbind(1:nrow(x),x))[missing.rows.indices,]

  return ( list (missing.matrix = missing.matrix,
                 numMissing = numMissing,
                 missing.rows.indices = missing.rows.indices,
                 missing.cols.indices = missing.cols.indices,
                 x.missing = x.missing) )
}
