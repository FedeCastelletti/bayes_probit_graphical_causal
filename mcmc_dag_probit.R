library(pcalg)
library(gRbase)
library(truncnorm)

source("move_dag_probit.r")
source("marg_like_dag.r")
source("posterior_sigma_dag_probit.r")
source("causal_effect_probit.r")


mcmc_dag_probit = function(y, X, TT, burn, a = NULL, g = NULL, w = NULL, causal = TRUE){
  
  ###########
  ## INPUT ##
  ###########
  
  # y :    (n,1) vector of observations from the binary response variable
  # X :    (n,p) matrix of observations from the p covariates
  # TT :   number of MCMC iterations
  # burn : burn-in period
  
  # a,g : hyperparamters of the DAG-Wishart prior
  # w :   prior probability of edge inclusion
  
  # causal : logical value (if TRUE, BMA causal effects are also computed; if FALSE, they are not)
  
  ############
  ## OUTPUT ##
  ############
  
  # A_chain :     matrix with TT columns, each representing the (vectorized) (q,q) adjacency matrix of one visited DAG
  # Sigma_chain : matrix with TT columns, each representing the (vectorized) (q,q) covariance matrix of one visited DAG
  # Causal_hat :  (n,q) matrix with BMA estimates of subject-specific causal effects (one row for each subject, one column for each intervened node)
  
  
  n = nrow(X)
  p = ncol(X)
  q = p + 1
  
  # Set DAG-Wishart hyperparameters (if not specified)
  
  if(is.null(a)){a = q}
  
  if(is.null(g)){g = 1/n}
  
  if(is.null(w)){w = 0.5}
  
  # Store space for parameters
  
  A_chain = matrix(0, nrow = q*q, ncol = TT)
  # matrix collecting the adjacecncy matrices of the DAGs
  # (each column is the by-column-vectorized adjacency matrix of the accepted DAG)
  
  Sigma_chain = matrix(0, nrow = q*q, ncol = TT)
  # matrix collecting the posterior draws from Sigma
  # (each column is the by-column-vectorized covariance matrix Sigma)
  
  L_chain = matrix(0, nrow = q*q, ncol = TT)
  D_chain = matrix(0, nrow = q*q, ncol = TT)
  
  z_chain = matrix(0, nrow = n, ncol = TT)
  # matrix collecting the posterior draws from the latent z
  
  Gamma_chain = rep(0, TT)
  # matrix collecting posterior samples for the cutoff gamma
  
  gamma_0 = 0
  
  inf_lim = rep(-Inf, n)
  sup_lim = rep(Inf, n)
  
  inf_lim[y == 1] = gamma_0
  sup_lim[y == 0] = gamma_0
  
  
  if(causal == TRUE){
    
    Causal_hat = matrix(0, n, q)
    # matrix collecting BMA causal effect estimates
    
  }else{
    
    Causal_hat = NULL
    
  }
  
  ## Inits values
  
  A_0 = matrix(0, q, q); colnames(A_0) = rownames(A_0) = 1:q; A_chain[,1] = A_0
  
  z_0 = rtruncnorm(n, a = inf_lim, b = sup_lim, mean = 0, sd = 1); z_chain[,1] = z_0
  
  Chol_Sigma_0 = posterior_sigma(cbind(z_0,X), A_0, g, a)
  
  Sigma_0 = Chol_Sigma_0$Sigma_post; Sigma_chain[,1] = Sigma_0
  
  L_0    = -Chol_Sigma_0$L_post[-1,1]
  D_1_0  = Chol_Sigma_0$D_post[1,1]
  
  A = A_0
  z = z_0
  
  gamma = gamma_0
  
  
  ## MCMC iterations
  
  cat("MCMC sampling")
  pb = utils::txtProgressBar(min = 2, max = TT, style = 3)
  
  for(t in 1:TT){
    
    ## update of DAG D
    
    A_move = move(A = A, q = q)
    
    A_star = A_move$A_new
    
    type.operator = A_move$type.operator
    
    nodes_star = A_move$nodes
    
    
    # Distinguish 3 cases:
    
    if(type.operator == 1){
      
      # (1) Insert a directed edge
      
      logprior = log(w/(1-w))
      
      marg_star = marg_like_j(j = nodes_star[2], dag = A_star, X = cbind(z,X), a = a, g = g, n = n)
      marg      = marg_like_j(j = nodes_star[2], dag = A, X = cbind(z,X), a = a, g = g, n = n)
      
    }else{
      
      if(type.operator == 2){
        
        # (2) Delete a directed edge
        
        logprior = log((1-w)/w)
        
        marg_star = marg_like_j(j = nodes_star[2], dag = A_star, X = cbind(z,X), a = a, g = g, n = n)
        marg      = marg_like_j(j = nodes_star[2], dag = A, X = cbind(z,X), a = a, g = g, n = n)
        
      }else{
        
        # (3) Reverse a directed edge
        
        logprior = log(1)
        
        marg_star = marg_like_j(j = nodes_star[1], dag = A_star, X = cbind(z,X), a = a, g = g, n = n) +
          marg_like_j(j = nodes_star[2], dag = A_star, X = cbind(z,X), a = a, g = g, n = n)
        
        marg = marg_like_j(j = nodes_star[1], dag = A, X = cbind(z,X), a = a, g = g, n = n) +
          marg_like_j(j = nodes_star[2], dag = A, X = cbind(z,X), a = a, g = g, n = n)
        
      }
      
    }
    
    # acceptance ratio
    
    ratio_D = min(0, marg_star - marg + logprior)
    
    # accept move
    
    if(log(runif(1)) < ratio_D){
      A = A_star
    }
    
    ## Sample Sigma
    
    Chol_Sigma = posterior_sigma(Y = cbind(z,X), A, g, a)
    Sigma_post = Chol_Sigma$Sigma_post
    L_post     = Chol_Sigma$L_post
    D_post     = Chol_Sigma$D_post
    
    ## Sample and update the latent
    
    Beta = -L_post[-1,1]
    
    z = rtruncnorm(n, a = inf_lim, b = sup_lim, mean = X%*%Beta, sd = sqrt(D_post[1,1]))
    
    # Update the cut-off
    
    g_prop = rnorm(1, gamma, 0.5)
    
    inf_lim_g = inf_lim
    inf_lim_g[y == 1] = g_prop
    sup_lim_g = sup_lim
    sup_lim_g[y == 0] = g_prop
    
    r_g = sum(log(pnorm(sup_lim_g, X%*%Beta, sd = sqrt(D_post[1,1])) - pnorm(inf_lim_g, X%*%Beta, sd = sqrt(D_post[1,1])))) -
      sum(log(pnorm(sup_lim, X%*%Beta, sd = sqrt(D_post[1,1])) - pnorm(inf_lim, X%*%Beta, sd = sqrt(D_post[1,1])))) +
      dnorm(g, gamma, 0.5, log = TRUE) - dnorm(gamma, g, 0.5, log = TRUE)
    
    # acceptance ratio
    
    ratio_g = min(0, r_g)
    
    # accept move
    
    if(log(runif(1)) < ratio_g){
      gamma   = g_prop
      inf_lim = inf_lim_g
      sup_lim = sup_lim_g
    }
    
    # store chain values
    
    A_chain[,t]     = A
    Sigma_chain[,t] = Sigma_post
    Gamma_chain[t]  = gamma
    
    
    ## compute causal effects (if causal = TRUE)
    
    if(causal == TRUE){
      
      if(t > burn){
        
        out_t = cbind(0, sapply(FUN = causal_y, X = 2:q, Sigma_D = Sigma_post, A = A, X_mat = X, g = gamma))
        
        Causal_hat = Causal_hat + out_t/(TT - burn)
        
      }
      
    }
    
    utils::setTxtProgressBar(pb, t)
    close(pb)
    
    
  }
  
  return(list(A_chain     = A_chain, 
              Sigma_chain = Sigma_chain,
              Causal_hat  = Causal_hat
  ))
}
