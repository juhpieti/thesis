//
functions {
  // function for log-likelihood for left-censored beta regression
  real lc_beta_density(real y, real mu, real rho, real a, real eps, int N) {
    real ll; //log-likelihood
    real y_scaled; //scaled y
    
    // resctrict mu to open interval (0,1)
    real mu_shifted;
    mu_shifted = fmin(fmax(mu,eps), 1-eps);
    
    // go separately through cases y=0, y=1, y \in (0,1)
    if (y==0) {
      ll = log(fmax(1e-20, beta_cdf(a/(a+1) | mu_shifted*rho, (1-mu_shifted)*rho))); //avoid log(0)
    } else if (y == 1) { // beta distribution cannot take y=1, shift using [y*(N-1)+0.5]/2 as proposed in the literature
      y_scaled = ((y*(N-1)+0.5)/N + a)/(1+a);
      ll = beta_lpdf(y_scaled | mu_shifted*rho, (1-mu_shifted)*rho) - log(a+1);
    } else {
      y_scaled = (y+a)/(a+1);
      ll = beta_lpdf(y_scaled | mu_shifted*rho, (1-mu_shifted)*rho) - log(a+1);
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
  
  // hierarchical priors for beta
  vector[n_var/2] mu_beta_1; // mean for first order terms
  vector<upper=0>[n_var/2] mu_beta_2; // mean for second order terms
  real mu_alpha; // mean for intercepts

  vector<lower=0>[n_var/2] s_beta_1; // standard deviation for first order terms
  vector<lower=0>[n_var/2] s_beta_2; // standard deviation for second order terms
  real<lower=0> s_alpha; // standard deviation for intercepts
  
  vector[J] alpha; // intercept terms
  vector<lower=0>[J] rho; // scale parameters for beta distribution
  
  // latent factors & species loadings
  matrix[n_obs_grid,n_f] Z_spat; // N(0,1) to create n_f spatially correlated latent factors
  matrix[N,n_f] Z; // n_f latent factors for each sampling location (non-spatial)
  matrix[n_f,J] Lambda; // species loadings
  
  // spatial latent factors
  vector<lower=0>[n_f] l; //length-scale parameters
}

transformed parameters{
  matrix[n_obs_grid,n_f] phi; // n_f spatially correlated latent factors for n_obs_grid (number of observed grid cells) locations
  
  for (k in 1:n_f) {
    matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L; // cholesky decomposition
    
    K = gp_exponential_cov(s_array,0.5,l[k]); // reduce magnitude to s = 0.5 for total variance of 1 after adding non-spatial Z with 0.5 variance as well.
    K = K + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability
    L = cholesky_decompose(K);
    
    phi[,k] = L*Z_spat[,k]; //follows GP(0,K)
  }
  
  matrix[N,J] Mu; // gather the mean parameters for beta distribution in a matrix
  for (n in 1:N) {
    Mu[n] = inv_logit(to_row_vector(alpha) + X[n]*append_row(beta_1,beta_2) + (Z[n] + P[n]*phi) * Lambda); // intercept + fixed effects + random effects
  }
}

model {
  // priors for coefficients
  for (v in 1:n_var/2) {
    beta_1[v] ~ normal(mu_beta_1[v],s_beta_1[v]);
    beta_2[v] ~ normal(mu_beta_2[v],s_beta_2[v]);
  }
  
  alpha ~ normal(mu_alpha, s_alpha);
  
  // hyperpriors
  mu_beta_1 ~ normal(0,1);
  mu_beta_2 ~ normal(0,1);
  mu_alpha ~ normal(0,2);
  
  //s_beta_1 ~ cauchy(0,1);
  //s_beta_2 ~ cauchy(0,1);
  
  s_beta_1 ~ cauchy(0,0.5);
  s_beta_2 ~ cauchy(0,0.5);
  s_alpha ~ cauchy(0,1);
  
  //s_beta_1 ~ student_t(4,0,1);
  //s_beta_2 ~ student_t(4,0,1);
  //s_beta_1 ~ exponential(1);
  //s_beta_2 ~ exponential(1);

  //alpha ~ normal(0,sqrt(2)); // intercepts
  rho ~ cauchy(0,sqrt(10)); // precision parameter
  
  // factor loadings
  to_vector(Z_spat) ~ normal(0,1); // for creating phi ~ GP(0,K)
  to_vector(Z) ~ normal(0,0.5); // reduce variance from 1 to 0.5 for total variance of 1 (after spatial RE added)
  to_vector(Lambda) ~ normal(0,1);
  
  // length-scale parameters
  l ~ cauchy(8,sqrt(50));
  //l ~ student_t(2,8,sqrt(50));
  
  // likelihood terms
  for (i in 1:N) {
    for (j in 1:J) {
      target += lc_beta_density(Y[i,j],Mu[i,j],rho[j],a,eps,N);
    }
  }
}

generated quantities {
  // log-likelihood for LOO-calculations
  vector[N] log_lik;

  for (i in 1:N) {
    real ll = 0;
    for (j in 1:J) {
      ll += lc_beta_density(Y[i,j],Mu[i,j],rho[j],a,eps,N);
    }
    log_lik[i] = ll;
  }
}