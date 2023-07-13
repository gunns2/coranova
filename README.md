# coranova

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package uses sample correlations to formally compare two or more polygenic scores in one or more populations.  \
Two implementations are available: parametric and non-parametric. We recommend the parametric implementation for continuous outcomes, and non-parametric for binary outcomes, or other non-normally distributed traits. 

Paper forthcoming.


## Installation

You can install the development version of coranova from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gunns2/coranova")
library(coranova)
```
## Guidelines



### Data Preparation

To compare polygenic scores across multiple population samples, create a list of data frames where each data frame contains the data from a distinct population sample. Each data frame should contain a column with the outcome variable and computed polygenic scores to be compared. The names of the columns should be shared across population sample data frames.

### Example



AFR:
```
  pheno       pgs1        pgs2        pgs3
  -1.4606535  -1.1363873  0.7988800   0.05945586
  -0.1498208  0.9218904   -1.8164597  -1.95811704
  3.5076853   0.2239545   1.5862969   -0.12250429
```

EUR:
```
  pheno       pgs1        pgs2        pgs3
  -3.249809   0.66336192  -2.470082   -0.5937651
  -2.754595   -0.72927645 2.431823    -0.1291377
  -1.698634   0.07946734  -1.129588   -0.9772898
```

To compare all three PGS across the two populations:
``` r
perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"))
```

To compare pgs1 and pgs2 in African population: \
This command will give confidence interval for difference in correlation between outcome and the two scores in the African sample.
``` r
perform_coranova_parametric(list(afr), "pheno", c("pgs1", "pgs2"))
```
## Warnings

Both PGS and outcome should be adjusted for principal components prior to analysis.

In simulations we find at least 1000 bootstrap and 1000 permutations are necessary to control type 1 error with binary outcomes. We have not tested other non-normal outcome distributions. 

## References
Bilker WB, Brensinger C, Gur RC. A Two Factor ANOVA-like Test for Correlated Correlations: CORANOVA. Multivariate Behav Res. 2004 Oct 1;39(4):565-94. doi: 10.1207/s15327906mbr3904_1. PMID: 26745459.

Olkin I, Finn JD. Testing correlated correlations. Psychological Bulletin. 1990 Sep;108(2):330.

## Contact
Sophia Gunn [gunns2@bu.edu](gunns2@bu.edu)
