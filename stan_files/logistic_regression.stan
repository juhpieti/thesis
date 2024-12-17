//
data {
  int<lower=0> N; // 
  int<lower=0> n_var; // number of variables
  array[N] int<lower=0,upper=1> y;
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
  y ~ bernoulli_logit(alpha + X*beta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha + X[n]*beta);
  }
}

