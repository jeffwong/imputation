#This file is deprecated

#' Robust SVD Imputation
#'
#' Imputation using a robust implementation of the SVD
#' Formulate the SVD as a biconvex minimization problem, and solve
#' it by alternating least squares.  When computing least squares,
#' use the least trimmed squares implementation, a robust
#' implementation that trims extreme values out
#' @param x a data frame or matrix where each row represents a different record
#' @param k the rank-k approximation to use for x
#' @param alpha the alpha passed into least trimmed squares
#' @param max.iters the number of times to alternate least trimmed squares
#' @param verbose if TRUE print status updates
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

#' CV for robustSVDImpute
#'
#' Cross Validation for robust SVD Imputation
#' Artificially erase some data and run robustSVDImpute multiple times,
#' varying k from 1 to k.max.  For each k, compute the RMSE on the subset of x
#' for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param k.max the largest rank used to approximate x
#' @param parallel runs each run for k = 1 to k = k.max in parallel.  Requires
#'   a parallel backend to be registered
cv.robustSVDImpute = function(x, k.max=floor(ncol(x)/2), parallel = T) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  if (parallel) {
    rmse = foreach (i=1:k.max, .combine = c, .packages = c('imputation')) %dopar% {
      x.imputed = robustSVDImpute(x.train, i, verbose=F)$x
      error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
      sqrt(mean(error^2))
    }
  }
  else {
    rmse = sapply(1:k.max, function(i) {
      x.imputed = robustSVDImpute(x.train, i, verbose=F)$x
      error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
      sqrt(mean(error^2))
    })
  }
  list(k = which.min(rmse), rmse = rmse[which.min(rmse)],
       k.full = 1:k.max, rmse.full = rmse)
}
