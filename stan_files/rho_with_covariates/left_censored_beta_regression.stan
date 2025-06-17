//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  real<lower=0> a; // zero-point as fixed (e.g. 1)
}

parameters {
  // coefficients for beta distribution mean
  vector[n_var/2] beta_1; // first order coefficients
  vector<upper=0>[n_var/2] beta_2;
  
  // coefficients for modeling variance of Beta
  vector[n_var/2] beta_rho_1;
  vector[n_var/2] beta_rho_2;
  
  real alpha;
  real alpha_rho;
  //real<lower=0> rho; // scale parameter for beta distribution
}

transformed parameters{
  vector[N] mu;
  mu = inv_logit(alpha + X*append_row(beta_1,beta_2));
  
  vector[N] logit_rho;
  logit_rho = alpha_rho + X*append_row(beta_rho_1,beta_rho_2);
}

model {
  // priors for coefficients
  beta_1 ~ normal(0,sqrt(2));
  beta_2 ~ normal(0,sqrt(2));
  beta_rho_1 ~ normal(0,sqrt(2));
  beta_rho_2 ~ normal(0,sqrt(1));

  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  alpha_rho ~ normal(0,sqrt(4));
  //rho ~ cauchy(0,sqrt(10));
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08); // restrict to open (0,1) interval
    
    // calculate the rho
    real rho_i;
    rho_i = 1 + exp(fmin(logit_rho[i],10));
    
    if (y[i] == 0) {
      target += beta_lcdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i);
      // target += beta_proportion_lcdf(a/(a+1) | mu_i, phi);
    } else if (y[i] ==  1) {
      // substract a small value
      (y[i]-1e-08+a)/(a+1) ~ beta(mu_i*rho_i, (1-mu_i)*rho_i);
      //target += log(1/(1+a))
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho_i, (1-mu_i)*rho_i); // scale to [0,1] interval
      //target += log(1/1+a)
    }
  }
}

generated quantities {
  // return beta_2 in constrained space
  //vector[n_var/2] beta_2;
  //beta_2 = -exp(logneg_beta_2);
  
  vector[N] log_lik;

  for (i in 1:N) {
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
    // calculate the rho
    real rho_i;
    rho_i = 1 + exp(fmin(logit_rho[i],10));
    
    if (y[i] == 0) {
      log_lik[i] = beta_lcdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i);
    } else if (y[i]==1) {
      log_lik[i] = beta_lpdf((y[i]-1e-08+a)/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i) - log(1+a);
    } else {
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i) - log(1+a); // scale to [0,1] interval
    }
  }
}
