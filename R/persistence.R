#' used for panel data in the wide format
persistenceImpute = function(x, verbose=T) {

    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) return (x)
    missing.matrix = prelim$missing.matrix
    missing.rows.indices = prelim$missing.rows.indices

    x[missing.rows.indices,] = sapply(missing.rows.indices, function(i) {
        bad.indices = which(missing.matrix[i,])
        good.indices = which(!missing.matrix[i,])
        x[i, bad.indices] = sapply(bad.indices, function(j) {
            #get nearest data point in the future
            neighbor = which.min(good.indices > j)
            if(length(neighbor) > 0) return (x[i,good.indices[neighbor]])
            else {
                neighbor = which.max(good.indices < j)
                if(length(neighbor) > 0) return (x[i,good.indices[neighbor]])
                    else {
                        warning("No good data to persist on")
                        return (NA)
                    }
            }
        })
        x[i,]
    })

    return ( list (
        x=x,
        missing.matrix=missing.matrix
    ))
}

cv.persistenceImpute = function(x) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.imputed = persistenceImpute(x.train, verbose=F, ...)$x
  error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
  rmse = sqrt(mean(error^2))
  
  list(imputation = x.imputed, rmse = rmse)
}
