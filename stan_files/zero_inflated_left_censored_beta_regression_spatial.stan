//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  int<lower=0> n_obs_grid; // number of observed grid cells
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  matrix[n_obs_grid,2] s; // coordinates of observed grid cells
  real<lower=0> a; // zero-point as fixed (e.g. 1)
  matrix[N,n_obs_grid] P; // binary matrix (nxm) telling to which grid cell each data point belongs 
}

parameters {
  // regression coefficients
  vector[n_var] beta;
  vector[n_var] beta_pi; // coefficients for probability of suitability
  real alpha;
  real alpha_pi; // intercept for probability of suitability
  
  // scale parameter for Beta distribution
  real<lower=0> rho; // scale parameter for beta distribution
  
  // GP Random Effect Parameters
  real<lower=0> l_mu; // length-scale parameter
  real<lower=0> s2_mu; // magnitude of the covariance function

  real<lower=0> l_pi; // length-scale parameter
  real<lower=0> s2_pi; // magnitude of the covariance function
  
  vector[n_obs_grid] z; // N(0,1) for creating the random effects
}

transformed parameters {
  vector[n_obs_grid] phi_mu; // random effects
  vector[n_obs_grid] phi_pi;
  {
    matrix[n_obs_grid,n_obs_grid] K_mu; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] K_pi;
    matrix[n_obs_grid,n_obs_grid] L_mu; // Cholesky of covariance matrix
    matrix[n_obs_grid,n_obs_grid] L_pi;
    
    // construct K
    for (i in 1:n_obs_grid) {
      for (j in i:n_obs_grid) {
        K_mu[i,j] = s2_mu*exp(-(distance(s[i],s[j]))/l_mu);
        K_mu[j,i] = s2_mu*exp(-(distance(s[i],s[j]))/l_mu);
        
        K_pi[i,j] = s2_pi*exp(-(distance(s[i],s[j]))/l_pi);
        K_pi[j,i] = s2_pi*exp(-(distance(s[i],s[j]))/l_pi);
      }
    }
    
    // generate the random effects
    L_mu = cholesky_decompose(K_mu);
    phi_mu = L_mu*z; // follows N(0,K)
    
    L_pi = cholesky_decompose(K_pi);
    phi_pi = L_pi*z;
  }
}

model {
  // priors for coefficients
  beta ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var))); // ORIGINAL
  beta_pi ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var)));
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  alpha_pi ~ normal(0,sqrt(2));
  rho ~ inv_gamma(1,1);
  
  z ~ normal(0,1);
  
  s2_mu ~ inv_gamma(0.25,0.25);
  s2_pi ~ inv_gamma(0.25,0.25);
  l_mu ~ inv_gamma(1,1);
  l_pi ~ inv_gamma(1,1);
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter for beta distribution
    real mu_i;
    mu_i = inv_logit(alpha + X[i,] * beta + P[i,]*phi_mu);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08); // there was problems of getting exact 0/1 with beta distribution calls
    // calculate the probability of suitability
    
    real pi_i;
    pi_i = inv_logit(alpha_pi + X[i,] * beta_pi + P[i,]*phi_pi);
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
    mu_i = inv_logit(alpha + X[i,] * beta + P[i,]*phi_mu);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = inv_logit(alpha_pi + X[i,] * beta_pi + P[i,]*phi_pi);
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

