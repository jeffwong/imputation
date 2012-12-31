cross.fn.default = function(x, i, j) {
  alpha = 0.5
  combo = c(alpha, 1 - alpha) %*% c(mean(x[i,]), mean(x[,j]))
  combo[[1]]
}

crossImpute = function(x, cross.fn = cross.fn.default) {
  prelim = impute.prelim(x)
  missing.matrix = prelim$missing.matrix
  missing.pairs = cbind(row(x)[missing.matrix], col(x)[missing.matrix])
  imputed = apply(missing.pairs, 1, function(i,j) {
    cross.fn(x,i,j)
  })
}
