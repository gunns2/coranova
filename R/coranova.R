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
#' \dontrun{mat <- matrix(c(1, 0.5, 0.6, 0.5, 1, 0.2, 0.6, 0.3, 1), ncol = 3)}
#' \dontrun{datA <- as.data.frame(MASS::mvrnorm(n = 100,  rep(0, 3), mat))}
#' \dontrun{datB <- as.data.frame(MASS:mvrnorm(n = 100,  rep(0, 3), mat))}
#' \dontrun{a <- list(cor(datA), cor(datB))}
#' \dontrun{perform_coranova(a, "V1", c("V2", "V3"), "parametric")}
#'
perform_coranova_parametric <- function(dat_list, outcome, measures){
  if(typeof(dat_list[[1]]) == "double"){
    dat_list <- list(dat_list)
  }
  cormat_list <- lapply(dat_list, cor)
  n_list <- lapply(dat_list, nrow)
  R <- populate_R(cormat_list, outcome, measures)
  V <- populate_V(cormat_list, outcome, measures, n_list)
  g <- length(dat_list)
  m <- length(measures)
  p <- length(measures) + 1
  if(g == 1){
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)
    if(m == 2){#if only two scores, can return difference and se for difference
      diff <- t(CW%*%R)
      se <- (1/solve(CW%*%V%*%t(CW)))^0.5
      return(list(pW = pW, diff = diff, se =  se,  LCB = diff - 1.96*se, UCB = diff + 1.96*se))
    }else{#if more than two scores, just return pval
      return(list(pW = pW))
    }
  }else if(m == 1){
    CB <- generate_linear_contrasts(g, m, "between")
    SB <- t(CB%*%R)%*%solve(CB%*%V%*%t(CB))%*%(CB%*%R)
    pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)
    if(g == 2){
      se <- (1/solve(CB%*%V%*%t(CB)))^0.5
      diff <- t(CB%*%R)
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
#' @examples \dontrun{TO ADD}
perform_coranova_nonparametric <- function(dat_list, outcome, measures, B, n_perm, test){
  perms <- c(rep(0, n_perm))
  nmeasures <- length(measures)
  method <- "bootstrap"
  print(nmeasures)
  if(test == "within"){
    stat <- "SW"
    for(i in 1:n_perm){
      #perform vars switching, within
      dat_list1 <- lapply(dat_list, shuffled_df, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list1, outcome, measures, method, B, stat)
    }
  }else if(test == "between"){
    test <- "SB"
    for(i in 1:n_perm){
      #perform group switching, between
      dat_list1 <- shuffle_groups(dat_list)
      perms[i] <- coranova_perm_helper(dat_list1, outcome, measures, method, B, stat)
    }
  }else if(test == "int"){
    test <- "SI"
    for(i in 1:n_perm){
      #first, group switching
      dat_list1 <- shuffle_groups(dat_list)
      #first, var switching
      dat_list2 <- lapply(dat_list1, shuffled_df, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list2, "V1", measures, method = method, B = B, stat = stat)
    }
  }
  chistat <- coranova_perm_helper(dat_list, "V1", measures, method = method, B = B, stat = stat)
  return(mean(as.numeric(chistat) <  perms))
}

