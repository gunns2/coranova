# coranova

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package can be used to formally compare two or more polygenic scores in one or more populations using correlation with outcome variable to compare the scores. 
Two implementations are available: parametric and non-parametric. We recommend the parametric implementation for continuous outcomes, and non-parametric for binary outcomes, or other non-normally distributed traits. 



## Installation

You can install the development version of coranova from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gunns2/coranova")
```
## Guidelines

In simulations we find at least 1000 bootstrap and 1000 permutations are necessary to control type 1 error with binary outcomes.

### Data Preparation

To compare polygenic scores across multiple population samples, create a list of data frames where each data frame contains the data from a distinct population sample. Each data frame should contain a column with the outcome variable and computed polygenic scores to be compared. The names of the columns should be shared across population sample data frames.


## Example


``` r
library(coranova)

head(afr)
head(eur)
perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"))

perform_coranova_parametric(list(afr), "pheno", c("pgs1", "pgs2"))
```

## References
Bilker WB, Brensinger C, Gur RC. A Two Factor ANOVA-like Test for Correlated Correlations: CORANOVA. Multivariate Behav Res. 2004 Oct 1;39(4):565-94. doi: 10.1207/s15327906mbr3904_1. PMID: 26745459.

Olkin I, Finn JD. Testing correlated correlations. Psychological Bulletin. 1990 Sep;108(2):330.

## Contact
Sophia Gunn [gunns2@bu.edu](gunns2@bu.edu)
