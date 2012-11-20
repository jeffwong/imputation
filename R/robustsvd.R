robustSVDImpute = function(x, k, alpha = 1/2, max.iters = 10, verbose = T ) {

    prelim = impute.prelim(x)
    if (prelim$numMissing == 0) return (x)
    missing.matrix = prelim$missing.matrix

    n = nrow(x)
    U = matrix(rnorm(n*k), n, k)

    for (i in 1:max.iters) {
        if(verbose) print(paste("Iteration:", i))
        DV.T = apply(x, 2, function(j) {
            ltsReg(U, j, intercept=F, alpha=alpha)$coefficients
        })
        VD = t(DV.T)
        U = t(apply(x, 1, function(i) {
            ltsReg(VD, i, intercept=F, alpha=alpha)$coefficients
        }))
    }
    
    UDVT = (U %*% DV.T)
    x[missing.matrix] = UDVT[missing.matrix]
    return ( list (
        x=x,
        missing.matrix=missing.matrix
    ))
}

cv.robustSVDImpute = function(x, k.max=floor(ncol(x)/2)) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  rmse = sapply(1:k.max, function(i) {
    x.imputed = robustSVDImpute(x.train, i, verbose=F)$x
    error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
    sqrt(mean(error^2))
  })
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse)
}
