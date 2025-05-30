# Factor augmented hurdle Poisson inverse model
Factor augmented hurdle Poisson inverse model (FAHP) can conduct supervised dimensionality reduction on high dimensional count data containing a large number of zeros.

# Installation
```r
install.packages("devtools")  
devtools::install_github("Wangzhijingwzj/fahp")  
library(fahp) 
```

# Usage
```r
fahp(zipdata, y, q = 2, maxit = 300, constraint = 5, para = 0)
```
* zipdata: Observation matrix, commonly characterized by zero inflation.
* y: A response vector with dimension equal to the number of rows in zipdata.
* q: The number of latent factors. Default as 2.
* maxit: Maximum number of iterations within optim function, defaults to 300.
* constraint: Constraint constant. Default as 5.
* para: The number of cores used for parallel processing. Defaults to 0, indicating that parallel processing is not performed.


# Value
The estimator results.
* FA: Factor score matrix.
* mu1: The intercept vector for hurdle part.
* mu2: The intercept vector for positive count part.
* t: Number of iterations.
* obj: The loglikelihhod results (omitting the constant), which can be used for model selection.


# Example
```r
n=50
p=20
x=matrix(rpois(n*p,lambda = 1),ncol = p)
y=rnorm(n)
fahp(x,y,q=1,maxit = 30,para = 10)
```
