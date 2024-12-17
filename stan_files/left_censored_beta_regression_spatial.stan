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
  vector[n_var] beta;
  // vector[n_var] beta_phi; // coefficients for phi parameter
  real alpha;
  // real alpha_phi; 
  real<lower=0> rho; // scale parameter for beta distribution
  //real alpha_phi;
  // real<lower=0> a; // zero-point parameter
  // vector<lower=0,upper=1>[N] V; // latent beta-distributed variables
  real<lower=0> l; // length-scale parameter
  real<lower=0> s2_cf; // magnitude of the covariance function
  vector[n_obs_grid] z; // N(0,1) for creating the random effects
}

transformed parameters {
  vector[n_obs_grid] phi; // random effects
  {
    matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L; // Cholesky of covariance matrix
    
    // construct K
    for (i in 1:n_obs_grid) {
      for (j in i:n_obs_grid) {
        K[i,j] = s2_cf*exp(-(distance(s[i],s[j]))/l);
        K[j,i] = s2_cf*exp(-(distance(s[i],s[j]))/l);
      }
    }
    // generate the random effects
    L = cholesky_decompose(K);
    phi = L*z; // follows N(0,K)
  }
}

model {
  // priors for coefficients
  beta ~ multi_normal(rep_vector(0,n_var),diag_matrix(rep_vector(sqrt(2),n_var))); // ORIGINAL
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  rho ~ inv_gamma(1,1);
  //alpha_phi ~ normal(0,1);
  //beta_phi ~ normal(0,1);
  
  z ~ normal(0,1);
  s2_cf ~ inv_gamma(0.25,0.25);
  l ~ inv_gamma(1,1);
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    mu_i = inv_logit(alpha + X[i,] * beta + P[i,]*phi);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    //real phi;
    //phi = exp(alpha_phi + X[i,] * beta_phi)+1e-08;
    if (y[i] == 0) {
      target += beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho);
      // target += beta_proportion_lcdf(a/(a+1) | mu_i, phi);
    } else if (y[i] ==  1) {
      // substract a small value
      (y[i]-1e-08+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho);
      //target += log(1/(1+a))
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho); // scale to [0,1] interval
      //target += log(1/1+a)
    }
  }
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    mu_i = inv_logit(alpha + X[i,] * beta + P[i,]*phi);
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
    //real phi;
    //phi = exp(alpha + X[i,] * beta_phi)+1e-08;
    if (y[i] == 0) {
      log_lik[i] = beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho);
      // log_lik[i] = beta_proportion_lcdf(a/(a+1) | mu_i, phi);
    } else if (y[i]==1) {
      log_lik[i] = beta_lpdf((y[i]-1e-08+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a);
    } else {
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a); // scale to [0,1] interval
      // log_lik[i] = beta_proportion_lpdf((y[i]+a)/(a+1) | mu_i, phi);
    }
  }
}
