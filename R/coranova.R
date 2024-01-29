#' Perform Coranova
#'
#' @param dat_list list of data frames, where each data frame refers to a separate population sample
#' @param outcome name of outcome variable (must be common across dataframes in dat_list)
#' @param measures names of measures to be compared
#'
#' @return results of coranova test, when group == 1 only performs within test, when nscores = 2 provides CI
#' @export
#'
#' @examples
#' perform_coranova_parametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"))
#'
perform_coranova_parametric <- function(dat_list, outcome, measures){
  if(typeof(dat_list[[1]]) == "double"){
    dat_list <- list(dat_list)
  }
  cormat_list <- lapply(dat_list, cor)
  n_list <- lapply(dat_list, nrow)
  R <- populate_mu(cormat_list, outcome, measures)
  V <- populate_sigma(cormat_list, outcome, measures, n_list)
  g <- length(dat_list)
  m <- length(measures)
  p <- length(measures) + 1
  if(g == 1){
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)
    if(m == 2){#if only two scores, can return difference and se for difference
      diff <- as.numeric(t(CW%*%R))
      se <- as.numeric((CW%*%V%*%t(CW))^0.5)
      return(list(pW = pW, diff = diff, se =  se,  LCB = diff - 1.96*se, UCB = diff + 1.96*se))
    }else{#if more than two scores, just return pval
      return(list(pW = pW))
    }
  }else if(m == 1){
    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)
    if(g == 2){
      se <- as.numeric((CB%*%V%*%t(CB))^0.5)
      diff <- as.numeric(t(CB%*%R))
      return(list(pB = pB, diff = diff, se = se, LCB = diff - 1.96*se, UCB = diff + 1.96*se))
    }else{
      return(list(pB = pB))
    }
  }else{
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)

    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)

    CI <-  generate_linear_contrasts(g, m, "interaction")
    SI <- t(CI%*%R)%*%solve(CI%*%V%*%t(CI))%*%(CI%*%R)
    pI <- pchisq(as.numeric(SI), (p-2)*(g-1),  lower.tail = F)
    return(list(pB = pB, pW = pW, pI = pI))
  }
}

#' Perform User-Specified Hypothesis Test
#'
#' @param dat_list list of data frames, where each data frame refers to a separate population sample
#' @param outcome name of outcome variable (must be common across dataframes in dat_list)
#' @param measures names of measures to be compared
#' @param contrast contrast matrix to generate hypothesis test
#' @param method to request parametric or bootstrap implementation of covariance matrix V
#' @param B number of bootstraps if method == "boot" is chosen
#'
#' @return results of coranova test, when group == 1 only performs within test, when nscores = 2 provides CI
#' @export
#'
#' @examples
#' #TO DO
perform_alt_test <- function (dat_list, outcome, measures, contrast, method, B)
{
  if (typeof(dat_list[[1]]) == "double") {
    dat_list <- list(dat_list)
  }
  cormat_list <- lapply(dat_list, cor)
  n_list <- lapply(dat_list, nrow)
  R <- populate_mu(cormat_list, outcome, measures)
  if (method == "parametric") {
    V <- populate_sigma(cormat_list, outcome, measures, n_list)
  }
  else if (method == "boot") {
    V <- bootstrap_sigma(dat_list, B, length(dat_list), outcome,
                     measures)
  }

  rank <- rankMatrix(contrast)
  S <- t(contrast %*% R) %*% solve(contrast %*% V %*% t(contrast)) %*% (contrast %*% R)

  p <- pchisq(as.numeric(S), rank, lower.tail = F)

  return(list(S = S, P = p))

}

#' Perform Coranova with bootstrapped V
#'
#' @param dat_list list of data frames, where each data frame refers to a separate population sample
#' @param outcome name of outcome variable
#' @param measures vector of measures to be compared
#' @param B number of bootstrap samples to run
#'
#' @return p-value of chosen test
#' @export
#'
#' @examples
#'  #perform_coranova_bootV(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"), 5, 5, "int")
#'
#'
perform_coranova_bootV <- function(dat_list, outcome, measures, B){
  if(typeof(dat_list[[1]]) == "double"){
    dat_list <- list(dat_list)
  }
  cormat_list <- lapply(dat_list, cor)
  n_list <- lapply(dat_list, nrow)
  R <- populate_mu(cormat_list, outcome, measures)
  V <- bootstrap_sigma(dat_list, B, length(dat_list), outcome, measures)
  g <- length(dat_list)
  m <- length(measures)
  p <- length(measures) + 1
  if(g == 1){
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)
    if(m == 2){#if only two scores, can return difference and se for difference
      diff <- as.numeric(t(CW%*%R))
      se <- as.numeric((CW%*%V%*%t(CW))^0.5)
      return(list(pW = pW, diff = diff, se =  se,  LCB = diff - 1.96*se, UCB = diff + 1.96*se))
    }else{#if more than two scores, just return pval
      return(list(pW = pW))
    }
  }else if(m == 1){
    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)
    if(g == 2){
      se <- as.numeric((CB%*%V%*%t(CB))^0.5)
      diff <- as.numeric(t(CB%*%R))
      return(list(pB = pB, diff = diff, se = se, LCB = diff - 1.96*se, UCB = diff + 1.96*se))
    }else{
      return(list(pB = pB))
    }
  }else{
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)

    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)

    CI <-  generate_linear_contrasts(g, m, "interaction")
    SI <- t(CI%*%R)%*%solve(CI%*%V%*%t(CI))%*%(CI%*%R)
    pI <- pchisq(as.numeric(SI), (p-2)*(g-1),  lower.tail = F)
    return(list(pB = pB, pW = pW, pI = pI))
  }
}


#' Perform Coranova with permutations
#'
#' @param dat_list list of data frames, where each data frame refers to a separate population sample
#' @param outcome name of outcome variable
#' @param measures vector of measures to be compared
#' @param B number of bootstrap samples to run
#' @param n_perm number of permutations to run
#' @param test type of test to be performed, options: "within", "between" or "int"
#'
#' @return p-value of chosen test
#' @export
#'
#' @examples
#'  perform_coranova_nonparametric(list(afr, eur), "pheno", c("pgs1", "pgs2", "pgs3"), 5, 5, "int")
#'
perform_coranova_nonparametric <- function(dat_list, outcome, measures, B, n_perm, test){
  perms <- c(rep(0, n_perm))
  nmeasures <- length(measures)
  if(test == "within"){
    stat <- "SW"
    for(i in 1:n_perm){
      #perform vars switching, within
      dat_list1 <- lapply(dat_list, shuffled_df, outcome = outcome, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list1, outcome, measures, B, stat)
    }
  }else if(test == "between"){
    stat <- "SB"
    for(i in 1:n_perm){
      #perform group switching, between
      dat_list1 <- shuffle_groups(dat_list)
      perms[i] <- coranova_perm_helper(dat_list1, outcome, measures, B, stat)
    }
  }else if(test == "int"){
    stat <- "SI"
    for(i in 1:n_perm){
      #first, group switching
      dat_list1 <- shuffle_groups(dat_list)
      #first, var switching
      dat_list2 <- lapply(dat_list1, shuffled_df,  outcome = outcome, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list2, outcome, measures, B = B, stat = stat)
    }
  }
  chistat <- coranova_perm_helper(dat_list, outcome, measures,  B = B, stat = stat)
  return(mean(as.numeric(chistat) <  perms))
}

