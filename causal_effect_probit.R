## Causal effect estimation

library(gRbase)

pa = function(set, object) {
  pa = which(object[,set] != 0)
  return(pa)
}

fa = function(set, object) {
  pa = which(object[,set] != 0)
  fa = c(set, pa)
  return(fa)
}

causal_y = function(Sigma_D, A, i, X_mat, g){
  
  # i : position of interveded node i in {2,...,q}
  # X is (n,p), where p = q - 1
  # The corresponding column in X is then (i-1)
  
  # return a (n,1) vector collecting the individual causal effects of do(X_j = x_j) on the binary response
  # defined as Pr(Y = 1 | do(X_j = x_j), Sigma_D, x_i)
  
  pa_D = pa(set = i, A)
  fa_D = fa(set = i, A)
  
  sigma_2_z = as.vector(Sigma_D[1,1] - Sigma_D[1,fa_D]%*%solve(Sigma_D[fa_D,fa_D])%*%Sigma_D[fa_D,1])
  
  gamma_D = Sigma_D[1,fa_D]%*%solve(Sigma_D[fa_D,fa_D])
  
  gamma_i  = gamma_D[1]
  gamma_pa = gamma_D[-1]
  
  if(length(pa_D) > 0){
    
    T_mat = solve(Sigma_D[pa_D,pa_D]) + (gamma_pa%*%t(gamma_pa))/sigma_2_z
    ss = sigma_2_z/(1 - t(gamma_pa)%*%solve(T_mat)%*%gamma_pa/sigma_2_z)
    
  }else{
    
    T_mat = 0
    ss = sigma_2_z
    
  }
  
  return(1 - pnorm(g, gamma_i*X_mat[,i-1], sqrt(ss)))
  
}
