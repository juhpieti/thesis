
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  int<lower=0> n_obs_grid; // number of observed spatial grid cells
  int<lower=0> J; // number of species
  int<lower=0> n_f; // number of latent factors
  //matrix<lower=0,upper=100>[N,J] Y; // species percentages
  array[N, J] int<lower=0, upper=100> Y; // species percentages
  //vector<lower=0>[N] y_sum; // total percent cover (sum(y))
  array[N] int<lower=0> y_sum; // total percent cover (sum(y))
  matrix[N,n_var] X; // covariate matrix
  matrix[n_obs_grid,2] s; // coordinates of observed grid cells
  matrix[N,n_obs_grid] P; // binary matrix (nxm) telling to which grid cell each data point belongs 
}

transformed data {
  real eps = 1e-12; // for y = 0
  array[n_obs_grid] vector[2] s_array; // kernel function takes input as array of vectors
  for(i in 1:n_obs_grid) {
    s_array[i] = to_vector(s[i]);
  }
}

parameters {
  // coefficients for Dirichlet distribution
  matrix[n_var/2,J] beta_1; // first order terms
  matrix<upper=0>[n_var/2,J] beta_2; // second order terms, restricted negative (bell-shaped curves based on ecological niche theory)
  vector[J] alpha; // intercept terms
  real<lower=0> scale_Dir; // scaling factor
  
  // latent factors & species loadings
  matrix[n_obs_grid,n_f] Z_spat; // N(0,1) to create n_f spatially correlated latent factors
  matrix[N,n_f] Z; // n_f latent factors for each sampling location (non-spatial)
  matrix[n_f,J] Lambda; // species loadings
  vector<lower=0>[n_f] l; //length-scale parameters

  // parameters for the total coverage model
  vector[n_var/2] beta_M_1; // first order terms
  vector<upper=0>[n_var/2] beta_M_2; // second order terms, restricted negative (bell-shaped curves based on ecological niche theory)
  real alpha_M; // intercept
  
  // spatial RE for the total coverage model
  vector[n_obs_grid] z_M; // N(0,1) to create random effects for total coverage model 
  real<lower=0> l_M; // length-scale
  real<lower=0> s2_M; // magnitude
}

transformed parameters{
  // Spatially correlated latent factors
  matrix[n_obs_grid,n_f] phi; // n_f spatially correlated latent factors for n_obs_grid (number of observed grid cells) locations
  
  for (k in 1:n_f) {
    matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L; // cholesky decomposition
    
    K = gp_exponential_cov(s_array,0.5,l[k]); // reduce magnitude to s = 0.5 for total variance of 1 after adding non-spatial Z with 0.5 variance as well.
    K = K + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability
    L = cholesky_decompose(K);
    
    phi[,k] = L*Z_spat[,k]; //follows GP(0,K)
  }
  
  // Parameters for the Dirichlet distribution
  matrix[N,J] Mu_Dir; // gather the mean parameters for beta distribution
  //array[N,J] real<lower=0,upper=1> Mu_Dir; // gather the mean parameters for beta distribution
  for (n in 1:N) {
    // Mu_Dir[n] = softmax(to_row_vector(alpha) + X[n]*append_row(beta_1,beta_2) + Z[n]*Lambda); // intercept + fixed effects + random effects
    Mu_Dir[n] = to_row_vector(softmax(to_vector(alpha) + to_vector(X[n]*append_row(beta_1,beta_2)) + to_vector( (Z[n] + P[n]*phi) * Lambda ))); // intercept + fixed effects + random effects
  }
  
  // Spatial RE for total coverage
  vector[n_obs_grid] phi_M;
  // Parameters for the model for total coverage
  vector[N] mu_M;
  
    {
    matrix[n_obs_grid,n_obs_grid] K_M; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L_M; // Cholesky of covariance matrix

    // construct K
    K_M = gp_exponential_cov(s_array,sqrt(s2_M),l_M);
    K_M = K_M + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability
    
    // generate the random effects
    L_M = cholesky_decompose(K_M);
    phi_M = L_M*z_M; // follows N(0,K)
    
    // mean parameter for Poisson
    mu_M = exp(alpha_M + X*append_row(beta_M_1,beta_M_2) + P*phi_M);
  }
}

model {
  // LIKELIHOODS
  for (i in 1:N) {
    // OBSERVATION MODEL
    Y[i] ~ dirichlet_multinomial(scale_Dir*to_vector(Mu_Dir[i]));
    // PROPORTIONS
    //Pi[i,] ~ dirichlet(scale_Dir*Mu_Dir[i,]);
    // TOTAL COVER
    y_sum[i] ~ poisson(mu_M[i]);
  }
  
  // PRIORS
  // Dirichlet component
  to_vector(beta_1) ~ normal(0,sqrt(2));
  to_vector(beta_2) ~ normal(0,sqrt(2));
  alpha ~ normal(0,10);
  scale_Dir ~ cauchy(0,sqrt(10)); // scaling/precision parameter for Dirichlet

  // Total Cover Component
  beta_M_1 ~ normal(0,sqrt(2));
  beta_M_2 ~ normal(0,sqrt(2));
  alpha_M ~ normal(0,10);
  
  // factor loadings (Dirichlet part)
  to_vector(Z_spat) ~ normal(0,1); // for creating phi ~ GP(0,K)
  to_vector(Z) ~ normal(0,0.5); // reduce variance from 1 to 0.5 for total variance of 1 (after spatial RE added)
  to_vector(Lambda) ~ normal(0,0.5);
  l ~ cauchy(8,sqrt(50)); // length-scale parameter
  
  // spatial RE (for total coverage)
  z_M ~ normal(0,1);
  l_M ~ cauchy(8,sqrt(50));
  s2_M ~ student_t(4,0,sqrt(0.1));
}

generated quantities {
  // log-likelihood for LOO-calculations
  vector[N] log_lik;

  for (i in 1:N) {
    log_lik[i] = poisson_lpmf(y_sum[i] | mu_M[i]) + dirichlet_multinomial_lpmf(Y[i] | scale_Dir*to_vector(Mu_Dir[i]));
  }
}
