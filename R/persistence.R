persistenceImpute = function(x, verbose=T) {
    missing.matrix = is.na(x)
    numMissing = sum(missing.matrix)
    if(verbose) {
        print(paste("imputing on", numMissing, "missing values with matrix size",
        nrow(x)*ncol(x), sep=" "))
    }
    if(numMissing == 0) {
        return (x)
    }

    missing.rows.indices = which(apply(missing.matrix, 1, function(i) {
        any(i)
    }))

    x[missing.rows.indices,] = sapply(missing.rows.indices, function(i) {
        bad.indices = which(missing.matrix[i,])
        good.indices = which(!missing.matrix[i,])
        x[i, bad.indices] = sapply(bad.indices, function(j) {
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
