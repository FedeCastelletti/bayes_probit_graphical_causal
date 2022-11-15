These R codes implement the Bayesian methodology of Castelletti & Consonni (2021, Bayesian Analysis) for structure learning and causal inference in probit graphical models.

Specifically:

mcmc_dag_probit.R : contains the main MCMC algorithm for posterior inference of DAG structures, parameters (covariance matrix) and causal effect estimation

move_dag_probit.R : performs one move from a DAG to an adjacent DAG (i.e. implements the proposal distribution over the space of DAGs satisfying the contraints of node 1 being the response)

marg_like_dag.R   : computes the (log)marginal likelihood of a DAG model relative to a node-component of the DAG

posterior_sigma_dag_probit: samples from a DAG-Wishart (posterior) distribution for Sigma; this is adapted for the contraint imposed by the probit function linking the continous covariates X to the binary response y (i.e. conditional variance of latent response (node 1) equal to 1)

causal_effect_probit.R : function to compute causal effects on response y following an interventions on a given Xj from an input covariance matrix Sigma

examples.R : implements mcmc_dag_probit.R and other auxiliary functions on simulated data
