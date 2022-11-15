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

m_j = function(j, dag, X, n, a, g){
  
  ## This function computes the marginal likelihood of a dag given the data X relative to a generic node j in {2,...,q}
  
  # n is the sample size
  # a,g are hyperparameters of the DAG-Wishart prior
  
  pa_j = pa(j, dag)
  
  y  = X[,j]
  XX = as.matrix(X[,pa_j])
  
  p_j   = length(pa_j)
  j_pos = q - p_j
  a_j_star = (a + q - 2*j_pos + 3)/2 - p_j/2 - 1
  
  if(length(pa_j) == 0){
    
    m = - 0.5*n*log(2*pi) + 
            lgamma(a_j_star + n/2) - lgamma(a_j_star) +
              a_j_star*log(0.5*g) - (a_j_star + n/2)*log(0.5*(g + sum(y^2)))
      
    
  } else{
    
    T_j = g*diag(1, p_j)
    
    T_j_hat = T_j + t(XX)%*%XX
    L_j_hat = solve(T_j_hat)%*%t(XX)%*%y
    
    m = - 0.5*n*log(2*pi) +
            0.5*log(det(T_j)) - 0.5*log(det(T_j_hat)) +
              lgamma(a_j_star + n/2) - lgamma(a_j_star) +
                a_j_star*log(0.5*g) - (a_j_star + n/2)*log(0.5*(g + sum(y^2) - t(L_j_hat)%*%T_j_hat%*%L_j_hat))
    
  }
  
  return(m)
  
}

m_1 = function(dag, X, n, a, g){
  
  ## This function computes the marginal likelihood of a dag given the data X relative to node 1 (the 0-1 response variable)
  
  # n is the sample size
  # a,g are hyperparameters of the DAG-Wishart prior
  
  pa_1 = pa(1, dag)
  
  y  = X[,1]
  XX = as.matrix(X[,pa_1])
  
  if(length(pa_1) == 0){
    
    m = -0.5*n*log(2*pi) -0.5*sum(y^2)
    
  } else{
    
    p_1 = length(pa_1)
    
    DD = t(XX)%*%XX
    
    V = DD + diag(g, p_1)
    
    b_hat = solve(V)%*%t(XX)%*%y
    e_hat = y - XX%*%b_hat
    
    m = - 0.5*n*log(2*pi) - 0.5*sum(y^2) +
      0.5*p_1*log(g) - 0.5*log(det(V)) + 0.5*t(b_hat)%*%V%*%b_hat
    
    
  }
  
  return(m)
  
}


marg_like_j = function(j, dag, X, a, g, n){
  
  ## General function to compute the marginal likelihood relative to a node j in the dag
  
  # n is the sample size
  # a,g are hyperparameters of the DAG-Wishart prior
  
  if(j != 1){
    
    out_m = m_j(j = j, dag = dag, X = X, n = n, a = a, g = g)
    
  }else{
    
    out_m = m_1(dag = dag, X = X, n = n, a = a, g = g)
    
  }
  
  return(out_m)
  
}
