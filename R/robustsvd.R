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
