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
  int<lower=0> J; // number of species
  int<lower=0> n_f; // number of latent factors
  matrix<lower=0,upper=1>[N,J] Y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  real<lower=0> a; // zero-point as fixed (e.g. 1)
}

transformed data {
  real eps = 1e-08; // small value for computational stability (avoid exact 0 and 1 in certain cases)
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
  matrix[N,n_f] Z; // n_f latent factors for each sampling location
  matrix[n_f,J] Lambda; // species loadings
  
  // latent factors & species loadings (suitability)
  matrix[N,n_f] Z_pi; // n_f latent factors for each sampling location
  matrix[n_f,J] Lambda_pi; // species loadings
}

transformed parameters{
  matrix[N,J] Mu; // gather the mean parameters for beta distribution
  for (n in 1:N) {
    Mu[n] = inv_logit(to_row_vector(alpha) + X[n]*append_row(beta_1,beta_2) + Z[n]*Lambda); // intercept + fixed effects + random effects
  }
  
  matrix[N,J] Pi; // gather the probability of suitability
  for (n in 1:N) {
    Pi[n] = inv_logit(to_row_vector(alpha_pi) + X[n]*append_row(beta_pi_1,beta_pi_2) + Z_pi[n]*Lambda_pi); // intercept + fixed effects + random effects
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
