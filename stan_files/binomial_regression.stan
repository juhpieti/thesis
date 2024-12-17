//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  array[N] int<lower=0,upper=100> y;
  matrix[N,n_var] X;
}

parameters {
  vector[n_var] beta;
  real alpha;
}

model {
  // priors for coefficients
  beta ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var))); // ORIGINAL
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  
  // likelihood
  y ~ binomial_logit(100, alpha + X*beta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | 100, alpha + X[n]*beta);
  }
}
