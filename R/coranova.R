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




#when group == 1 only performs within test, when nscores = 2 provides CI
perform_coranova <- function(dat_list, outcome, measures, method, B){
  if(typeof(dat_list[[1]]) == "double"){
    dat_list <- list(dat_list)
  }
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




perform_coranova_perm <- function(dat_list, outcome, measures, n, method, B, n_perm, test, stat){
  perms <- c(rep(0, n_perm))
  nmeasures <- length(measures)
  print(nmeasures)
  if(test == "within" & stat == "SW"){
    for(i in 1:n_perm){
      #perform vars switching, within
      dat_list1 <- lapply(dat_list, shuffled_df, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list1, n, outcome, measures, method, B, stat)
    }
  }else if(test == "between" & stat == "SB"){
    for(i in 1:n_perm){
      #perform group switching, between
      dat_list1 <- shuffle_groups(dat_list)
      perms[i] <- coranova_perm_helper(dat_list1, n, outcome, measures, method, B, stat)
    }
  }else if(test == "int" & stat == "SI"){
    for(i in 1:n_perm){
      #first, group switching
      dat_list1 <- shuffle_groups(dat_list)
      #first, var switching
      dat_list2 <- lapply(dat_list1, shuffled_df, measures = measures, nmeasures = nmeasures)
      perms[i] <- coranova_perm_helper(dat_list2, n, "V1", measures, method = method, B = B, stat = stat)
    }
  }
  chistat <- coranova_perm_helper(dat_list, n, "V1", measures, method = method, B = B, stat = stat)
  return(mean(as.numeric(chistat) <  perms))
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
coranova_perm_helper <- function(dat_list, n, outcome, measures, method, B, stat){
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
    #pB <- pchisq(as.numeric(SB), g-1,  lower.tail = F)
    return(SB)}
  if(stat == "SW"){
    CW <- generate_linear_contrasts(g, m, "within")
    SW <-  t(CW%*%R)%*%solve(CW%*%V%*%t(CW))%*%(CW%*%R)
    #pW <- pchisq(as.numeric(SW), p-2,  lower.tail = F)
    return(SW)}
  if(stat == "SI"){
    CI <-  generate_linear_contrasts(g, m, "interaction")
    SI <- t(CI%*%R)%*%solve(CI%*%V%*%t(CI))%*%(CI%*%R)
    #pI <- pchisq(as.numeric(SI), (p-2)*(g-1),  lower.tail = F)
    return(SI)}
}
