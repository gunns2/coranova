
#if two variables are the same must be listed first
get_cov_from_cor_mat <- function(cor_mat, vars1, vars2, n){
  val <- (0.5*(2*cor_mat[vars1[2], vars2[2]] - (cor_mat[vars1[1], vars1[2]]*cor_mat[vars2[1], vars2[2]]))*(1 - cor_mat[vars1[2], vars2[2]]^2 - cor_mat[vars1[1], vars1[2]]^2 - cor_mat[vars2[1], vars2[2]]^2 ) + cor_mat[vars1[2], vars2[2]]^3 )/n
  return(val)
}

generate_linear_contrasts <- function(n_groups, n_measures, type){
  g <- n_groups
  p <- n_measures + 1 #n measures + outcome variable
  if(type == "between"){
    C1 = kronecker(rep(1, g-1), t(rep(1, p -1 )))
    C2 = kronecker(diag(g-1), t(rep(-1,p-1)))
    C = cbind(C1, C2)
    return(C)
  }else if(type == "within"){
    C = kronecker(t(rep(1, g)), cbind(rep(1, p-2), -1*diag(p-2) ))
    return(C)
  }else if(type == "interaction"){
    C1 = kronecker(rep(1, g-1), cbind(rep(1, p-2), -1*diag(p-2)))
    C2 = kronecker(diag(g-1), cbind(rep(-1, p -2), diag(p -2)))
    C = cbind(C1, C2)
    return(C)
  }
}


#' Populate Sample Correlation Vector
#'
#' @param cormat_list list of correlation matrices, one for each sample
#' @param outcome name of outcome variable
#' @param measures vector of names of measure variables
#'
#' @return vector of correlations with measures and outcome variable
#' @export
#'
#' @examples
#' \dontrun{mat <- matrix(c(1, 0.5, 0.6, 0.5, 1, 0.2, 0.6, 0.3, 1), ncol = 3)}
#' \dontrun{datA <- as.data.frame(MASS::mvrnorm(n = 100,  rep(0, 3), mat))}
#' \dontrun{datB <- as.data.frame(MASS:mvrnorm(n = 100,  rep(0, 3), mat))}
#' \dontrun{a <- list(cor(datA), cor(datB))}
#' \dontrun{populate_R(a, "V1", c("V2", "V3"))}
populate_R <- function(cormat_list, outcome, measures){
  R <- lapply(cormat_list, getVals, outcome = outcome, measures = measures)
  return(unlist(R))
}

getVals <- function(cormat, outcome, measures){
  return(cormat[outcome, measures])
}

populate_V_one_group <- function(cor_mat, n, outcome, measures){
  V <- matrix(nrow = length(measures), ncol = length(measures))
  for (i in 1:length(measures)){
    for(j in  1:length(measures)){
      V[i,j] <- get_cov_from_cor_mat(cor_mat, c(outcome, measures[i]), c(outcome, measures[j]), n)
    }
  }
  return(V)
}

populate_V <- function(cormat_list, outcome, measures, n_list){
  stopifnot(length(cormat_list) == length(n_list))
  cov_list <- mapply(populate_V_one_group, cormat_list, n_list, MoreArgs = list(outcome = outcome, measures = measures), SIMPLIFY = F)
  return(bdiag(cov_list))
}

sample_dat <- function(dat){
  return(dat[sample(nrow(dat), replace = TRUE),])
}

bootstrap_V <- function(dat_list, B, num_pops, outcome, measures){
  stopifnot(num_pops == length(dat_list))
  num_scores <- length(measures)
  #step 3 repeat steps 1 & 2
  R_df <- matrix(nrow = B, ncol = num_scores*num_pops)

  for(b in 1:B){
    #step 1
    #Select a bootstrap sample from the original data set, using the following procedure: Randomly select a bootstrap sample, with replacement,
    #from group 1, consisting of n1 observations. Repeat this process for each group, selecting bootstrap samples from groups 2, ..., g
    #of sizes n2, ..., ng
    boot_list <- lapply(dat_list, sample_dat)
    #step 2
    #Compute the correlations in the vector R based on the bootstrap sample.
    cormat_boot_list <- lapply(boot_list, cor)
    R_df[b,]<- populate_R(cormat_boot_list, outcome, measures)
  }
  #step 4
  #Compute the covariance matrix R based on the B bootstrapped Rs, ˆ Boot V , which is a bootstrap estimate of V.
  return(round(cov(R_df), 5))
}

#permutation functions
shuffle_groups <- function(dat_list){
  dat_df <- bind_rows(dat_list, .id = "dat")
  dat_df1 <- dat_df %>% mutate(dat = sample(dat_df$dat, length(dat_df$dat), replace = F))
  return(split(dat_df1, f = dat_df1$dat, drop = T) %>% map(~select(.x, !c(`dat`))))
}

shuffle_vars <- function(row) {
  row[sample(length(row))]
}

shuffled_df <- function(dat, measures, nmeasures){
  shuffled <- apply(dat[, names(dat) %in% c(measures)], 1, shuffle_vars)
  dat1 <- cbind(dat$V1, as.data.frame(t(shuffled)))
  colnames(dat1) <- paste0("V", 1:c(nmeasures +1))
  return(dat1)
}

#coranova to be performed within permutation step to generate p-values, returns test
#statistic instead of p value
coranova_perm_helper <- function(dat_list,outcome, measures, method, B, stat){
  cormat_list <- lapply(dat_list, cor)
  n_list <- lapply(dat_list, nrow)
  R <- populate_R(cormat_list, outcome, measures)
  if(method == "parametric"){
    V <- populate_V(cormat_list, outcome, measures, n_list)
  }else if(method == "boot"){
    V <- bootstrap_V(dat_list, B, length(dat_list), outcome, measures)
  }
  g <- length(dat_list)
  m <- length(measures)
  p <- length(measures) + 1
  if(stat == "SB"){
    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    return(SB)}
  if(stat == "SW"){
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    return(SW)}
  if(stat == "SI"){
    CI <-  generate_linear_contrasts(g, m, "interaction")
    SI <- t(CI%*%R)%*%solve(CI%*%V%*%t(CI))%*%(CI%*%R)
    return(SI)}
}


