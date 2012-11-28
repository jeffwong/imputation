\name{gbmImpute}
\alias{gbmImpute}
\title{GBM Imputation}
\usage{
  gbmImpute(x, max.iters = 2, cv.fold = 2, n.trees = 100,
    verbose = T, ...)
}
\arguments{
  \item{x}{a data frame or matrix where each row is a
  different record}

  \item{max.iters}{number of times to iterate through the
  columns and impute each column with fitted values from a
  regression tree}

  \item{cv.fold}{number of folds that gbm should use
  internally for cross validation}

  \item{n.trees}{the number of trees used in gradient
  boosting machines}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation using Boosted Trees Fill each column by
  treating it as a regression problem.  For each column i,
  use boosted regression trees to predict i using all other
  columns except i.  If the predictor variables also
  contain missing data, the gbm function will itself use
  surrogate variables as substitutes for the predictors.
  This imputation function can handle both categorical and
  numeric data.
}

