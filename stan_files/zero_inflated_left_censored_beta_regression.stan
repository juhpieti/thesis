//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  real<lower=0> a; // zero-point as fixed (e.g. 1)
}

parameters {
  vector[n_var] beta;
  // vector[n_var] beta_phi; // coefficients for phi parameter
  vector[n_var] beta_pi; // coefficients for probability of suitability
  real alpha;
  // real alpha_phi; // intercept for phi parameter
  real alpha_pi; // intercept for probability of suitability
  real<lower=0> rho; // scale parameter for beta distribution
  // real<lower=0,upper=1> pi_i; // probability of suitability (this will have to be modeled by covariates)
  // real<lower=0> a; // zero-point parameter
}

model {
  // priors for coefficients
  beta ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var))); // ORIGINAL
  beta_pi ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var)));
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  alpha_pi ~ normal(0,sqrt(2));
  rho ~ inv_gamma(1,1);
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter for beta distribution
    real mu_i;
    mu_i = inv_logit(alpha + X[i,] * beta);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08); // there was problems of getting exact 0/1 with beta distribution calls
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = inv_logit(alpha_pi + X[i,] * beta_pi);
    pi_i = fmin(fmax(pi_i,1e-08), 1 - 1e-08);
    
    
    // likelihood terms
    if (y[i] == 0) {
      // logSumExp as in https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      target += log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
      // target += beta_proportion_lcdf(a/(a+1) | mu_i, phi);
    } else if (y[i]==1) {
      (y[i]-1e-08+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho);
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho); // scale to [0,1] interval
      // (y[i]+a)/(a+1) ~ beta_proportion(mu_i, phi);
      target += log(pi_i);
    }
  }
}


generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    mu_i = inv_logit(alpha + X[i,] * beta);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
    // calculate the probability of suitability
    
    real pi_i;
    pi_i = inv_logit(alpha_pi + X[i,] * beta_pi);
    pi_i = fmin(fmax(pi_i,1e-08), 1 - 1e-08);
    
    
    // likelihood terms
    if (y[i] == 0) {
      log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
      // log_lik[i] = beta_proportion_lcdf(a/(a+1) | mu_i, phi);
    } else if (y[i] == 1) {
      log_lik[i] = beta_lpdf((y[i]-1e-08+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a);
    } else {
      // I suspect the log of jacobian should be added! log(1/a+1)=log(1/2)=-0.69...!
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a); // scale to [0,1] interval
      // log_lik[i] = beta_proportion_lpdf((y[i]+a)/(a+1) | mu_i, phi);
      log_lik[i] += log(pi_i);
    }
  }
}


