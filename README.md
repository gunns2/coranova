# coranova

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package can be used to formally evaluate two ore more polygenic scores in one or more populations.  

## Installation

You can install the development version of coranova from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gunns2/coranova")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(coranova)
library(MASS)

mat <- matrix(c(1, 0.5, 0.6, 0.5, 1, 0.2, 0.6, 0.3, 1), ncol = 3)
datA <- as.data.frame(MASS::mvrnorm(n = 100,  rep(0, 3), mat))
datB <- as.data.frame(MASS:mvrnorm(n = 100,  rep(0, 3), mat))
perform_coranova(list(cor(datA), cor(datB)) , "V1", c("V2", "V3"), "parametric")
```

