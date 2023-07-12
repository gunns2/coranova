# coranova

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package can be used to formally compare two or more polygenic scores in one or more populations using correlation with outcome variable to compare the scores. 
Two implementations are available: parametric and non-parametric. We recommend the parametric implementation for continuous outcomes, and non-parametric for binary outcomes, or other non-normal traits. 

## Installation

You can install the development version of coranova from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gunns2/coranova")
```

## Example


``` r
library(coranova)
library(MASS)

mat <- matrix(c(1, 0.5, 0.6, 0.5, 1, 0.2, 0.6, 0.3, 1), ncol = 3)
datA <- as.data.frame(MASS::mvrnorm(n = 100,  rep(0, 3), mat))
datB <- as.data.frame(MASS:mvrnorm(n = 100,  rep(0, 3), mat))
perform_coranova(list(cor(datA), cor(datB)) , "V1", c("V2", "V3"), "parametric")
```

## References
Bilker WB, Brensinger C, Gur RC. A Two Factor ANOVA-like Test for Correlated Correlations: CORANOVA. Multivariate Behav Res. 2004 Oct 1;39(4):565-94. doi: 10.1207/s15327906mbr3904_1. PMID: 26745459.

Olkin I, Finn JD. Testing correlated correlations. Psychological Bulletin. 1990 Sep;108(2):330.

##Contact
Sophia Gunn [gunns2@bu.edu](gunns2@bu.edu)
