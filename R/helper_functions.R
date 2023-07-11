
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




