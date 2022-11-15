################################
## Examples on simulated data ##
################################

library(pcalg)

q = 20   # number of nodes
n = 1000 # sample size

# Generate true DAG and parameters

set.seed(1)

true_dag = randomDAG(q, prob = 3/(2*q-2))

A = t(as(true_dag, "matrix"))
A[A != 0] = 1

L = A*matrix(runif(q*q, 1, 2), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE); diag(L) = 1
D = diag(rep(1, q))

Sigma = solve(t(L))%*%D%*%solve(L)
mu    = c(rep(0, q))

# Generate the data

library(mvtnorm)

Y = rmvnorm(n, mu, Sigma)

m = colMeans(Y)
s = apply(X = Y, FUN = sd, MARGIN = 2)

Y = t(t(Y) - m)

head(Y)

# Example 1 : Compute marginal likelihood relative to one node (e.g. j = 1)

source("marg_like_dag.R")

marg_like_j(j = 1, dag = A, X = Y, a = q, g = 1, n = n)


# Example 2 : Perform local moves between DAGs

source("move_dag_probit.R")

move(A = A, q = q)


# Example 3 (sample from the posterior of Sigma)

source("posterior_sigma_dag_probit.R")

a = q
g = 1/n

sigma_post = posterior_sigma(Y, A, g, a)$Sigma_post

round(sigma_post, 2)


## Finally, create the 0-1 response variable by thresholding of Y[,1]

X = Y[,-1]
y = Y[,1]
y[y > 0] = 1
y[y < 0] = 0


## Data consists of y (binary response) and X (n,p) matrix with covariates


#################
## MCMC scheme ##
#################

## Fix number of MCMC iterations and burn-in period

TT   = 6000
burn = 1000


source("mcmc_dag_probit.R")

## Posterior inference on DAGs and parameters (Sigma)

t_0 = proc.time()
out = mcmc_dag_probit(y = y, X = X, TT = TT, burn = burn, causal = FALSE)
t_1 = proc.time() - t_0

## Posterior inference on DAGs and parameters (Sigma) + causal effect estimation

t_0 = proc.time()
out_causal = mcmc_dag_probit(y = y, X = X, TT = TT, burn = burn)
t_1 = proc.time() - t_0


#################################
## Compute posterior summaries ##
#################################

## Posterior probabilities of edge inclusion

probs = round(matrix(rowMeans(out_causal$A_chain[,(burn + 1):TT]), q, q), 4)
probs

## Median Probability DAG Model

A_hat = round(probs)
A_hat

## BMA estimate of the covariance matrix

Sigma_hat = round(matrix(rowMeans(out_causal$Sigma_chain[,(burn + 1):TT]), q, q), 2)
Sigma_hat

## BMA estimate of causal effects

round(out_causal$Causal_hat, 2)

