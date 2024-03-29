% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coranova.R
\name{perform_coranova_nonparametric}
\alias{perform_coranova_nonparametric}
\title{Perform Coranova with permutations}
\usage{
perform_coranova_nonparametric(dat_list, outcome, measures, B, n_perm, test)
}
\arguments{
\item{dat_list}{list of data frames, where each data frame refers to a separate population sample}

\item{outcome}{name of column containing outcome variable (must be common across dataframes in dat_list)}

\item{measures}{vector of column names of measures to be compared}

\item{B}{number of bootstrap samples to run}

\item{n_perm}{number of permutations to run}

\item{test}{type of test to be performed, options: "within", "between" or "int"}
}
\value{
p-value of chosen test
}
\description{
Perform Coranova with permutations
}
\examples{
 perform_coranova_nonparametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"), 5, 5, "int")

}
