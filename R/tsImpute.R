#Time series imputation using cubic splines
tsImpute = function(x, period) {
    n = nrow(x)
    apply(x, 2, function(j) {
        spline(1:n, j)$y
    })
}
