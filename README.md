imputation
==========

[imputation on CRAN](http://cran.r-project.org/web/packages/imputation/index.html)

Missing data imputation (also known as matrix completion) is an extremely difficult science that tries
to fill in missing values of a dataset with the best guess.  Recently, it was popularized by the Netflix Challenge,
where a matrix of Netflix users and their movie ratings were presented to the data science community to see
if algorithms could be developed to predict how a user would rate a certain movie that the user has not yet seen.

References:
* [Missing value estimation methods for DNA microarrays](http://bioinformatics.oxfordjournals.org/content/17/6/520.full.pdf).  Troyanskaya, et al.
* [A Singular Value Thresholding Algorithm for Matrix Completion](http://arxiv.org/pdf/0810.3286v1.pdf).  Cai, Candes, Shen.

##Imputation Algorithms Presented

* Mean Imputation
* k-Nearest Neighbors 
* SVD Imputation
* SVT Imputation
* Boosted Trees Imputation
* Locally weighted least squares

##Highlights

* meanImpute is a good way to start any missing data problem.  It's the fastest imputation technique and does reasonably well
* Sometimes, we want to identify missing values and impute them by fitting a line through its neighbors.  This can be done by taking a set of points {y_t, x_t} and regressing y_t on the index t.  Additionally, we can use a locally weighted least squares regression line to taylor the weights of the data points that are observed near the missing ones.  This is done in lmImpute
* gbmImpute is a technique to impute missing data when both categorical and numerical data is available.  It uses boosted decision trees, which requires lots of data in order to work well.  It has the advantage though of partitioning data, and then fitting different means to the partitions
* tsImpute is a technique to impute time series data.  There are three significant components to any time series problem: time, dimensions, and metrics.  The dimensions are categorical variables describing the data points, and metrics are the actual time series data.  tsImpute projects the time variable using [TimeProjection](https://github.com/jeffwong/TimeProjection), and then imputes the metrics using boosted trees again.  The time projections help to further segment the data points, for example identifying day vs night segments, weekday vs weekend segments, etc.
* kNN and SVD impute are classic imputation methods described in Troyanskaya.  The SVD finds a low rank k approximation to the data, which can be suitable for noisy data.  kNN is only good when the number of features is small  
* SVT is a recently popularized imputation algorithm that does very well with numeric data.  It is however the slowest algorithm presented here, requiring the computation of many SVDs.  SVTApproxImpute can be used as an estimate, simply computing the SVD once, thresholding the singular values at lambda, then multiplying the decomposition again to get the imputation

##Algorithm Design

Each function in this package includes the imputation algorithm as well as a cross validatiion algorithm.  The CV
algorithm artificially eliminates 1/3 of the data in a dataset, and runs the imputation function.  Using the completed
data, the RMSE is calculated on the portion of the data that was artificially removed only.  Different imputation
algorithms will perform differently on different datasets, so it is important to have these functions for comparison.

