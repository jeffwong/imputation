getSVDImpute = function(x) {
    k = 10
    num.iters = 10
    verbose = F
    SVDImpute(x, k = k, num.iters = num.iters, verbose=F)$x
}

getRobustSVDImpute = function(x) {
    k = 3
    alpha = 1
    max.iters = 3
    robustSVDImpute(x, k = k, alpha = alpha, max.iters = max.iters)$x
}

getSVTImpute = function(x) {
    SVTImpute(x, lambda = 0.1, verbose=F)$x
}

getKNNImpute = function(x) {
    kNNImpute(x, k = 3, verbose = F)$x
}

getgbmImpute = function(x) {
    gbmImpute(x, max.iters = 1, verbose = T)$x
}

getPersistenceImpute = function(x) {
    persistenceImpute(x, verbose = F)$x
}

getMeanImpute = function(x) {
    apply(x, 2, function(j) {
        bad.indices = which(is.na(j))
        j[bad.indices] = mean(j, na.rm = T)
        j
    })
}

#RANDOM DATA IMPUTATION

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getSVDImpute,
                                        numRow = 1000, numCol = 10,
                                        numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getRobustSVDImpute,
                                        numRow = 1000, numCol = 10,
                                        numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getSVTImpute,
                                        numRow = 1000, numCol = 10,
                                        numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getKNNImpute,
                                        numRow = 2000, numCol = 10,
                                        numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getgbmImpute,
                                        numRow = 1000, numCol = 10,
                                        numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.random(imputation.fn = getMeanImpute,
                                        numRow = 1000, numCol = 10,
                                        numMissing = 500))

#TIME SERIES IMPUTATION

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getSVDImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getRobustSVDImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getSVTImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getKNNImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getgbmImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getPersistenceImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))

set.seed(100)
system.time(imputation.benchmark.ts(imputation.fn = getMeanImpute,
                                    numTS = 1000, TSlength = 10,
                                    numMissing = 500))
