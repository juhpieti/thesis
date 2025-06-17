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

transformed data {
  real eps = 1e-08; // small value for computational stability (avoid exact 0 and 1 in certain cases)
}

parameters {
  // coefficients for beta distribution mean
  vector[n_var/2] beta_1; // first order coefficients
  vector<upper=0>[n_var/2] beta_2; //

  // coefficients for probability of suitability
  vector[n_var/2] beta_pi_1;
  vector<upper=0>[n_var/2] beta_pi_2;

  // intercepts
  real alpha;
  real alpha_pi;
  
  // scale parameter for Beta distribution
  real<lower=0> rho; // scale parameter for beta distribution
  
  // GP Random Effect Parameters
  real<lower=0> l_mu; // length-scale parameter
  real<lower=0> s2_mu; // magnitude of the covariance function

  real<lower=0> l_pi; // length-scale parameter
  real<lower=0> s2_pi; // magnitude of the covariance function
  
  vector[n_obs_grid] z; // N(0,1) for creating the random effects
  vector[n_obs_grid] z_pi;
  //vector[n_obs_grid] phi_mu; // spatial random effects
  //vector[n_obs_grid] phi_pi; //
}

transformed parameters {
  vector[n_obs_grid] phi_mu; // random effects
  vector[n_obs_grid] phi_pi;
  
  vector[N] mu;
  vector[N] prob_suit;

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
    
    K_mu = K_mu + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability
    K_pi = K_pi + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability

    // generate the random effects
    L_mu = cholesky_decompose(K_mu);
    phi_mu = L_mu*z; // follows N(0,K)

    L_pi = cholesky_decompose(K_pi);
    phi_pi = L_pi*z_pi;
    
    // mean for beta distribution
    mu = inv_logit(alpha + X*append_row(beta_1,beta_2) + P*phi_mu);
    // probability of suitability
    prob_suit = inv_logit(alpha_pi + X*append_row(beta_pi_1,beta_pi_2) + P*phi_pi);
  }
    
}

model {
  // construct K matrices
  //matrix[n_obs_grid,n_obs_grid] K_mu; // covariance matrix
  //matrix[n_obs_grid,n_obs_grid] K_pi;
  //matrix[n_obs_grid,n_obs_grid] L_mu; // Cholesky of covariance matrix
  //matrix[n_obs_grid,n_obs_grid] L_pi;

  // construct K
  // for (i in 1:n_obs_grid) {
  //   for (j in i:n_obs_grid) {
  //     K_mu[i,j] = s2_mu*exp(-(distance(s[i],s[j]))/l_mu);
  //     K_mu[j,i] = s2_mu*exp(-(distance(s[i],s[j]))/l_mu);
  // 
  //     K_pi[i,j] = s2_pi*exp(-(distance(s[i],s[j]))/l_pi);
  //     K_pi[j,i] = s2_pi*exp(-(distance(s[i],s[j]))/l_pi);
  //   }
  // }
  
  // K_mu = K_mu + diag_matrix(rep_vector(1e-08,n_obs_grid));
  // K_pi = K_pi + diag_matrix(rep_vector(1e-08,n_obs_grid));
  // 
  // L_mu = cholesky_decompose(K_mu);
  // L_pi = cholesky_decompose(K_pi);

  //phi_mu ~ multi_normal(rep_vector(0,n_obs_grid),K_mu); //GP(0,K)
  //phi_mu ~ multi_normal_cholesky(rep_vector(0,n_obs_grid),L_mu);
  //phi_pi ~ multi_normal(rep_vector(0,n_obs_grid),K_pi); 
  //phi_pi ~ multi_normal_cholesky(rep_vector(0,n_obs_grid),L_pi);

  // priors for coefficients
  beta_1 ~ normal(0,sqrt(1));
  beta_2 ~ normal(0,sqrt(1));
  beta_pi_1 ~ normal(0,sqrt(1));
  beta_pi_2 ~ normal(0,sqrt(1));

  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  alpha_pi ~ normal(-1,sqrt(0.25));
  
  rho ~ cauchy(0,sqrt(10));
  
  z ~ normal(0,1);
  z_pi ~ normal(0,1);

  // first try
  //s2_mu ~ cauchy(0,sqrt(0.01));
  //s2_pi ~ cauchy(0,sqrt(0.001));

  // second try
  //s2_mu ~ student_t(4,0,sqrt(0.1));
  //s2_pi ~ student_t(10,0,sqrt(0.1));
  
  // third try
  s2_mu ~ student_t(2,0,sqrt(0.1));
  s2_pi ~ student_t(4,0,sqrt(0.1));
  
  // fourth try
  //sqrt(s2_mu) ~ student_t(4,0,sqrt(0.1));
  //target += log(2*sqrt(s2_mu));
  //sqrt(s2_pi) ~ student_t(4,0,sqrt(0.1));
  //target += log(2*sqrt(s2_mu));
  
  l_mu ~ cauchy(8,sqrt(50));
  //l_pi ~ cauchy(8,sqrt(50));
  l_pi ~ student_t(2,8,sqrt(50));
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter for beta distribution
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,eps), 1 - eps); // there was problems of getting exact 0/1 with beta distribution calls
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = prob_suit[i];
    pi_i = fmin(fmax(pi_i,eps), 1 - eps);
    
    // likelihood terms
    if (y[i] == 0) {
      // logSumExp as in https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      //target += log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
      target += log_sum_exp(log(1-pi_i), log(pi_i) + log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho))));
      
    } else if (y[i]==1) { // stan declares beta on open interval (0,1)
      (y[i]-eps+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho);
      target += log(pi_i);
      
    // y \in (0,1)
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho); // scale to [0,1] interval
      target += log(pi_i);
    }
  }
}


generated quantities {
  // log-likelihood for LOO-calculations
  vector[N] log_lik;
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,eps), 1 - eps);
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = prob_suit[i];
    pi_i = fmin(fmax(pi_i,eps), 1 - eps);
    
    // likelihood terms
    if (y[i] == 0) {
      // log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
      log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho))));
    } else if (y[i] == 1) {
      log_lik[i] = beta_lpdf((y[i]-eps+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a) + log(pi_i);
    } else {
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a) + log(pi_i); // scale to [0,1] interval
    }
  }
}

