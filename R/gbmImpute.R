gbmImpute = function(x, max.iters = 2, cv.fold = 2, n.trees = 100, verbose=T, ...) {
  
  prelim = impute.prelim(x)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix
  missing.cols.indices = prelim$missing.cols.indices
  
  print(paste("Training over:", length(missing.cols.indices), "features"))
  for (i in 1:max.iters) {
    if (verbose) print (paste("Begin iteration: ", i))
    x[,missing.cols.indices] = sapply(missing.cols.indices, function(j) {
      if (verbose) print( paste("Imputing on feature: ", j))
      #response variable in gbm cannot contain missing data
      good.data = which(!missing.matrix[,j])
      gbm1 <- gbm(x[good.data,j] ~ .,
                  data = as.data.frame(x[good.data,-j]),
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
      data.predict = predict(gbm1, newdata = as.data.frame(x[-good.data,-j]), n.trees = best.iter)
      x[-good.data,j] = data.predict
      x[,j]
    })
  }
  
  return ( list (
    x=x,
    missing.matrix=missing.matrix
  ))
}

cv.gbmImpute = function(x, ...) {
  prelim = cv.impute.prelim(x)
  remove.indices = prelim$remove.indices
  x.train = prelim$x.train

  x.imputed = gbmImpute(x.train, verbose=F, ...)$x
  error = (x[remove.indices] - x.imputed[remove.indices]) / x[remove.indices]
  rmse = sqrt(mean(error^2))
  
  list(imputation = x.imputed, rmse = rmse)
}
