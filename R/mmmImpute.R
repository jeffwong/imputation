.Mode <- function(x) {
  ux <- unique(x)
  ux <- ux[!is.na(ux)]
  ux[which.max(tabulate(match(x, ux)))]
}

#' Mean-Mode-Median Imputation
#'
#' Fill missing values in a column with the mean, mode, or median of the column
#' depending on the data type
#' @param x a data frame or matrix where each row represents a different record
#' @examples
#'   x = matrix(rnorm(100),10,10)
#'   x.missing = x > 1
#'   x[x.missing] = NA
#'   mmmImpute(x)
#' @export
mmmImpute = function(x) {
  prelim = impute.prelim(x, byrow=F, keep.x.missing = F)
  if (prelim$numMissing == 0) return (x)
  missing.matrix = prelim$missing.matrix

  x.imputed = apply(x[,prelim$missing.cols.indices], 2, function(j) {
    bad.indices = which(is.na(j))
    datatype = data.class(j)
    if (datatype %in% c('factor', 'Date', 'character', 'logical')) {
      j[bad.indices] = .Mode(j)
    }
    else if (datatype == 'ordered') {
      #median of ordered factors.
      #When the sample has an even number of values, use the left value in the order
      #Arbitrarily chosen!
      j[bad.indices] = factor(levels(j)[floor(median(as.numeric(j)))], levels = levels(j), ordered=T)
    }
    else {
      j[bad.indices] = mean(j, na.rm=T)
    }
    return (j)
  })
  x[,prelim$missing.cols.indices] = x.imputed
  
  return (list(x = x,
               missing.matrix = missing.matrix))
}
