//
data {
  int<lower=0> N; // 
  int<lower=0> n_var; // number of variables
  array[N] int<lower=0,upper=1> y;
  matrix[N,n_var] X;
}

parameters {
  vector[n_var/2] beta_1;
  vector<upper=0>[n_var/2] beta_2;
  real alpha;
}

model {
  // priors for coefficients
  beta_1 ~ multi_normal(rep_vector(0,n_var/2),diag_matrix(rep_vector(sqrt(2),n_var/2))); // ORIGINAL
  beta_2 ~ multi_normal(rep_vector(0,n_var/2),diag_matrix(rep_vector(sqrt(2),n_var/2))); // ORIGINAL
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  
  // likelihood
  y ~ bernoulli_logit(alpha + X[,1:n_var/2]*beta_1 + X[,(1+n_var/2):n_var]*beta_2);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | alpha + X[n,1:n_var/2]*beta_1 + X[n,(1+n_var/2):n_var]*beta_2);
  }
}

