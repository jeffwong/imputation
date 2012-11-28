\name{persistenceImpute}
\alias{persistenceImpute}
\title{Persistence Imputation}
\usage{
  persistenceImpute(x, verbose = T)
}
\arguments{
  \item{x}{a data frame or matrix where each row contains
  one period's worth of measurements and columns represent
  times of that period}

  \item{verbose}{if TRUE print status updates}
}
\description{
  Imputation by filling missing values with the nearest
  known value in time. Each row of x is one period, missing
  values are filled by searching within the row and finding
  the nearest column which has data. Used for panel data in
  the wide format
}

