#' SVT Imputation
#'
#' Imputation using Singular Value Thresholding a la Cai, Candes, Shen.
#' @param x a data frame or matrix with size n1 x n2 where each row represents a different record
#' @param lambda the penalty on the singular values
#' @param stepsize optional.  If not provided, uses 1.2 * (n1 * n2) / (number of missing elements)
#' @param threshold convergence threshold
#' @param max.iters maximum number of iterations.  Note that each iteration will require computing
#'   an SVD
#' @param verbose if TRUE print status updates
#' @references A Singular Value Thresholding Algorithm for Matrix Completion. Cai, Candes, Shen. 
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   SVTImpute(x, 3)
#' @export
SVTImpute = function(x, lambda, stepsize, threshold = 1e-3, max.iters = 10, verbose=F) {
  prelim = impute.prelim(x, byrow=F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix

  x.zeroimpute = x
  x.zeroimpute[missing.matrix] = 0
  x.zeroimpute.norm = norm(x.zeroimpute, "F")

  if (missing(stepsize)) stepsize = min(1.2 * length(x) / sum(missing.matrix), 1.9)
  k = ceiling(lambda / ( stepsize * norm(x.zeroimpute, "F")))

  Y = k*stepsize*x.zeroimpute

  for (iter in 1:max.iters) {
    if (verbose) print(paste("Begin iteration", iter))
    y.svd = svd(Y)
    lambda.indices = which(y.svd$d < lambda)
    if(length(lambda.indices) > 0) {
      d.augmented = c(y.svd$d[-lambda.indices] - lambda, 
                      rep(0, length(lambda.indices)))
    } else {
      d.augmented = y.svd$d - lambda
    }
    X.k = (y.svd$u %*% diag(d.augmented) %*% t(y.svd$v))
    X.k.temp = X.k; X.k.temp[missing.matrix] = 0
    if (norm(X.k.temp - x.zeroimpute, "F") / x.zeroimpute.norm < threshold) {
      if (verbose) print(paste("Converging on iteration", iter))
      break
    }
    Y[!missing.matrix] = Y[!missing.matrix] + stepsize*(x[!missing.matrix] - X.k[!missing.matrix])
    Y[missing.matrix] = 0
  }

  return( list(
    x=X.k,
    missing.matrix = missing.matrix
  ))
}

#' CV for SVTImpute
#'
#' Cross Validation for SVT Imputation
#' Artificially erase some data and run SVTImpute multiple times,
#' varying lambda using lambda.range.  For each lambda, compute the
#' RMSE on the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param lambda.range a vector of penalty terms to use in the CV
#' @param parallel runs each run for lambda in lambda.range in parallel.  Requires
#'   a parallel backend to be registered
#' @param ... extra parameters to pass to SVTImpute
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   cv.SVTImpute(x)
#' @export
cv.SVTImpute = function(x, lambda.range = seq(0,1,length.out=101), parallel = F, ...) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  if (parallel) {
    if (!require(foreach)) stop("R package foreach is required for parallel execution, as well
                                 as a registered parallel backend")
    rmse = foreach (i=lambda.range, .combine = c, .packages = c('imputation')) %dopar% {
      x.imputed = SVTImpute(x.train, i, verbose=F, ...)$x
      error = (x.imputed[remove.indices] - x[remove.indices])
      nerror = error / x[remove.indices]
      list(nrmse = sqrt(mean(nerror^2)), rmse = sqrt(mean(error^2))) 
    }
  }
  else {
    rmse = sapply(lambda.range, function(i) {
      x.imputed = SVTImpute(x.train, i, verbose=F, ...)$x
      error = (x.imputed[remove.indices] - x[remove.indices])
      nerror = error / x[remove.indices]
      list(nrmse = sqrt(mean(nerror^2)), rmse = sqrt(mean(error^2)))
    })
  }
  nrmse = unlist(rmse[1,]); rmse = unlist(rmse[2,])
  list(lambda = lambda.range[which.min(rmse)], rmse = rmse[which.min(rmse)],
       lambda.full = lambda.range, rmse.full = rmse, nrmse.full = nrmse)
}
