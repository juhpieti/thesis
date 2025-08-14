//
functions {
  // function for log-likelihood for left-censored beta regression
  real zi_lc_beta_density(real y, real mu, real pi,real rho, real a, real eps, int N) {
    real ll; //log-likelihood
    real y_scaled; //scaled y
    
    // resctrict mu to open interval (0,1)
    real mu_shifted;
    mu_shifted = fmin(fmax(mu,eps), 1-eps);
    
    real pi_shifted;
    pi_shifted = fmin(fmax(pi,eps), 1 - eps);
    
    // go separately through cases y=0, y=1, y \in (0,1)
    if (y==0) {
      // logSumExp as in https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      ll = log_sum_exp(log(1-pi_shifted), log(pi_shifted) + log(fmax(1e-20, beta_cdf(a/(a+1) | mu_shifted*rho, (1-mu_shifted)*rho)))); //avoid log(0)
    } else if (y == 1) { // beta distribution cannot take y=1, shift using [y*(N-1)+0.5]/2 as proposed in the literature
      y_scaled = ((y*(N-1)+0.5)/N + a)/(1+a);
      ll = beta_lpdf(y_scaled | mu_shifted*rho, (1-mu_shifted)*rho) - log(a+1) + log(pi_shifted);
    } else {
      y_scaled = (y+a)/(a+1);
      ll = beta_lpdf(y_scaled | mu_shifted*rho, (1-mu_shifted)*rho) - log(a+1) + log(pi_shifted);
    }
    return ll;
  }
}

data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  int<lower=0> n_obs_grid; // number of observed grid cells
  int<lower=0> J; // number of species
  int<lower=0> n_f; // number of latent factors
  matrix<lower=0,upper=1>[N,J] Y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  matrix[n_obs_grid,2] s; // coordinates of observed grid cells
  real<lower=0> a; // zero-point as fixed (e.g. 1)
  matrix[N,n_obs_grid] P; // binary matrix (nxm) telling to which grid cell each data point belongs 
}

transformed data {
  real eps = 1e-08; // small value for computational stability (avoid exact 0 and 1 in certain cases)
  array[n_obs_grid] vector[2] s_array; //kernel function takes input as array of vectors
  for(i in 1:n_obs_grid) {
    s_array[i] = to_vector(s[i]);
  }
}

parameters {
  // coefficients for beta distribution mean
  matrix[n_var/2,J] beta_1; // first order terms
  matrix<upper=0>[n_var/2,J] beta_2; // second order terms, restricted negative (bell-shaped curves based on ecological niche theory)
  
  // coefficients for suitability
  matrix[n_var/2,J] beta_pi_1;
  matrix[n_var/2,J] beta_pi_2;
  
  vector[J] alpha; // intercept terms
  vector[J] alpha_pi; // intercepts for suitability
  vector<lower=0>[J] rho; // scale parameters for beta distribution
  
  // latent factors & species loadings (mean of beta)
  matrix[n_obs_grid,n_f] Z; // N(0,1) to create n_f latent factors for mean of beta
  matrix[n_f,J] Lambda; // species loadings
  
  // latent factors & species loadings (suitability)
  matrix[n_obs_grid,n_f] Z_pi; // N(0,1) to create n_f latent factors for location suitability
  matrix[n_f,J] Lambda_pi; // species loadings
  
  // lengthscale parameters for spatial latent factors
  vector<lower=0>[n_f] l; //length-scale for beta mean
  vector<lower=0>[n_f] l_pi; //length-scale for suitability  
}

transformed parameters{
  matrix[n_obs_grid,n_f] phi; // n_f spatially correlated latent factors for n_obs_grid (number of observed grid cells) locations
  matrix[n_obs_grid,n_f] phi_pi; // same for suitability of location
  
  for (k in 1:n_f) {
    matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L; // cholesky decomposition
    
    matrix[n_obs_grid,n_obs_grid] K_pi; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L_pi; // cholesky decomposition
    
    K = gp_exponential_cov(s_array,1,l[k]); // use s = 1 for variance of 1 as in non-spatial latent factor model
    L = cholesky_decompose(K);
    
    K_pi = gp_exponential_cov(s_array,1,l_pi[k]); // use s = 1 for variance of 1 as in non-spatial latent factor model
    L_pi = cholesky_decompose(K_pi);
    
    phi[,k] = L*Z[,k]; //follows GP(0,K)
    phi_pi[,k] = L_pi*Z_pi[,k];
  }
  
  matrix[N,J] Mu; // gather the mean parameters for beta distribution in a matrix
  for (n in 1:N) {
    Mu[n] = inv_logit(to_row_vector(alpha) + X[n]*append_row(beta_1,beta_2) + P[n]*phi*Lambda); // intercept + fixed effects + random effects
  }
  
  matrix[N,J] Pi;
  for (n in 1:N) {
    Pi[n] = inv_logit(to_row_vector(alpha_pi) + X[n]*append_row(beta_pi_1,beta_pi_2) + P[n]*phi_pi*Lambda_pi); // intercept + fixed effects + random effects
  }
}

model {
  // priors for coefficients
  to_vector(beta_1) ~ normal(0,sqrt(1));
  to_vector(beta_2) ~ normal(0,sqrt(1));
  
  to_vector(beta_pi_1) ~ normal(0,sqrt(1));
  to_vector(beta_pi_2) ~ normal(0,sqrt(1));

  alpha ~ normal(0,sqrt(2)); // intercepts
  alpha_pi ~ normal(0,sqrt(2));
  rho ~ cauchy(0,sqrt(10)); // precision parameter
  
  // factor loadings
  to_vector(Z) ~ normal(0,1);
  to_vector(Lambda) ~ normal(0,1);
  
  to_vector(Z_pi) ~ normal(0,1);
  to_vector(Lambda_pi) ~ normal(0,1);
  
  // lengthscale parameters
  l ~ cauchy(8,sqrt(50));
  l_pi ~ student_t(2,8,sqrt(50));
  
  // likelihood terms
  for (i in 1:N) {
    for (j in 1:J) {
      target += zi_lc_beta_density(Y[i,j],Mu[i,j],Pi[i,j],rho[j],a,eps,N);
    }
  }
}

generated quantities {
  // log-likelihood for LOO-calculations
  vector[N] log_lik;

  for (i in 1:N) {
    real ll = 0;
    for (j in 1:J) {
      ll += zi_lc_beta_density(Y[i,j],Mu[i,j],Pi[i,j],rho[j],a,eps,N);
    }
    log_lik[i] = ll;
  }
}
