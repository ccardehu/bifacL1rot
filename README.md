
# bifacL1rot

<!-- badges: start -->
<!-- badges: end -->

L1 rotation of Exploratory Generalized Bi-Factor Models

## Installation

You can install the development version of `bifacL1rot` from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("ccardehu/bifacL1rot")
```

## Example

This is a basic example:

``` r
set.seed(1234)
A <- matrix(0, 20, 4) # A 20 x 4 Factor loading matrix
A[,1] = runif(20, 1, 2) # Bi-factor structure, first factor
A[,2] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
A[,3] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
A[,4] = runif(20, 0.5, 1) * rbinom(20, 1, 0.25)
# Create rotation matrix (via random matrix):
Ah =  array(rnorm(length(A), sd = .5), dim = dim(A))
Tr = eigen(t(Ah) %*% Ah)$vectors
Ah = (A)%*%(Tr)
res = bifacL1rot::bifactorL1(Ah)
Arot = res$B
Phi = res$Phi
```

