#' GBM Imputation
#'
#' Imputation using Boosted Trees
#' Fill each column by treating it as a regression problem.  For each
#' column i, use boosted regression trees to predict i using all other
#' columns except i.  If the predictor variables also contain missing data,
#' the gbm function will itself use surrogate variables as substitutes for the predictors.
#' This imputation function can handle both categorical and numeric data.
#' @param x a data frame or matrix where each row is a different record
#' @param max.iters number of times to iterate through the columns and
#'   impute each column with fitted values from a regression tree
#' @param cv.fold number of folds that gbm should use internally for cross validation
#' @param n.trees the number of trees used in gradient boosting machines
#' @param verbose if TRUE print status updates
#' @param ... additional params passed to gbm
#' @examples
#'   x = matrix(rnorm(10000),1000,10)
#'   x.missing = x > 2
#'   x[x.missing] = NA
#'   gbmImpute(x)
#' @export
gbmImpute = function(x, max.iters = 2, cv.fold = 2, n.trees = 100, verbose = T, ...) {
  if (nrow(x) < 1000) warning("Tree based imputation works best with larger data (> 1000 obs)")
  prelim = impute.prelim(x, byrow = F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  missing.cols.indices = prelim$missing.cols.indices
 
  if (verbose) print(paste("Training over:", length(missing.cols.indices), "features"))
  for (i in 1:max.iters) {
    if (verbose) print (paste("Begin iteration: ", i))
    x[,missing.cols.indices] = sapply(missing.cols.indices, function(j) {
      if (verbose) print(paste("Imputing on feature: ", j))
      #response variable in gbm cannot contain missing data
      if (i == 1) good.data = which(!missing.matrix[,j])
      else good.data = 1:nrow(x)
      bad.data = which(missing.matrix[,j])
      data = x[good.data,-j]
      if (!is.data.frame(data)) data = data.frame(data = data)
      gbm1 <- gbm(x[good.data,j] ~ .,
                  data = data,
                  var.monotone = rep(0, ncol(x)-1), # -1: monotone decrease,
                  # +1: monotone increase,
                  #  0: no monotone restrictions
                  distribution="gaussian",     # bernoulli, adaboost, gaussian,
                  # poisson, coxph, and quantile available
                  n.trees=n.trees,                # number of trees
                  shrinkage=0.005,             # shrinkage or learning rate,
                  # 0.001 to 0.1 usually work
                  interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
                  bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
                  train.fraction = 2/3,        # fraction of data for training,
                  n.minobsinnode = 10,         # minimum total weight needed in each node
                  cv.folds = cv.fold,                # do 5-fold cross-validation
                  keep.data=TRUE,              # keep a copy of the dataset with the object
                  verbose=verbose,
                  ...)                # print out progress
      best.iter <- gbm.perf(gbm1,method="test", plot.it = F)
      newdata = x[bad.data,-j]
      if (!is.data.frame(newdata)) newdata = data.frame(data = newdata)
      data.predict = predict(gbm1, newdata = newdata, n.trees = best.iter)
      x[bad.data,j] = data.predict
      x[,j]
    })
  }
  
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}

#' CV for gbmImpute
#'
#' Cross Validation for GBM Imputation
#' Artificially erase some data and run gbmImpute.  Compute the RMSE
#' on the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param ... extra parameters to be passed to gbmImpute
#' @examples
#'   x = matrix(rnorm(10000),1000,10)
#'   x.missing = x > 2
#'   x[x.missing] = NA
#'   cv.gbmImpute(x)
#' @export
cv.gbmImpute = function(x, ...) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.imputed = gbmImpute(x.train, verbose=F, ...)$x
  error = (x[remove.indices] - x.imputed[remove.indices])
  nerror = error / x[remove.indices]
  rmse = sqrt(mean(error^2))
  nrmse = sqrt(mean(nerror^2))
  
  list(imputation = x.imputed, rmse = rmse, nrmse = nrmse)
}
