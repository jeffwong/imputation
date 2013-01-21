#' Time Series Imputation
#'
#' Time Series Imputation using Boosted Trees
#' Fill each column by treating it as a regression problem.  For each
#' column i, use boosted regression trees to predict i using all other
#' columns except i.  If the predictor variables also contain missing data,
#' the gbm function will itself use surrogate variables as substitutes for the predictors.
#' This imputation function can handle both categorical and numeric data.
#' @param time a vector of dates or datetime objects
#' @param dimension a data frame of exogenous predictor variables
#' @param metric a matrix where each column represents a time series
#' @param max.iters number of times to iterate through the columns and
#'   impute each column with fitted values from a regression tree
#' @param cv.fold number of folds that gbm should use internally for cross validation
#' @param n.trees the number of trees used in gradient boosting machines
#' @param verbose if TRUE print status updates
#' @export
tsImpute = function(time, dimension, metric, max.iters = 2, cv.fold = 2,
                    n.trees = 100, verbose = T, ...) {
  time.projection = projectDate(time)
  fixed = cbind(time.projection, dimension)
  prelim = impute.prelim(metric, byrow = F)
  if (prelim$numMissing == 0) return (metric)
  missing.matrix = prelim$missing.matrix
  missing.cols.indices = prelim$missing.cols.indices
 
  if (verbose) print(paste("Training over:", length(missing.cols.indices), "features"))
  for (i in 1:max.iters) {
    if (verbose) print (paste("Begin iteration: ", i))
    metric[,missing.cols.indices] = sapply(missing.cols.indices, function(j) {
      if (verbose) print(paste("Imputing on feature: ", j))
      #response variable in gbm cannot contain missing data
      if (i == 1) good.data = which(!missing.matrix[,j])
      else good.data = 1:nrow(metric)
      bad.data = which(missing.matrix[,j])
      gbm1 <- gbm(metric[good.data,j] ~ .,
                  data = as.data.frame(cbind(fixed[good.data,], metric[good.data,-j])),
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
                  train.fraction = 0.5,        # fraction of data for training,
                  # first train.fraction*N used for training
                  n.minobsinnode = 10,         # minimum total weight needed in each node
                  cv.folds = cv.fold,                # do 5-fold cross-validation
                  keep.data=TRUE,              # keep a copy of the dataset with the object
                  verbose=T)                # print out progress
      best.iter <- gbm.perf(gbm1,method="OOB", plot.it = F)
      data.predict = predict(gbm1,
                             newdata = as.data.frame(cbind(fixed[bad.data,], metric[bad.data,-j])),
                             n.trees = best.iter)
      metric[bad.data,j] = data.predict
      metric[,j]
    })
  }
  
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}

#' CV for tsImpute
#'
#' Cross Validation for Time Series Imputation
#' Artificially erase some data and run gbmImpute.  Compute the RMSE
#' on the subset of x for which data was artificially erased.
#' @param x a data frame or matrix where each row represents a different record
#' @param ... extra parameters to be passed to gbmImpute
#' @export
cv.gbmImpute = function(date, dimension, metric, ...) {
  prelim = cv.impute.prelim(metric)
  remove.indices = prelim$remove.indices
  metric.train = prelim$x.train
  date.train = date[-remove.indices]
  dimension.train = dimension[-remove.indices,]

  metric.imputed = tsImpute(date.train, dimension.train, metric.train, verbose=F, ...)$x
  error = (metric[remove.indices] - metric.imputed[remove.indices]) / metric[remove.indices]
  rmse = sqrt(mean(error^2))
  
  list(imputation = x.imputed, rmse = rmse)
}
