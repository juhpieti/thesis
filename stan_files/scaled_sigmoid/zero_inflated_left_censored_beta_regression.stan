//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  real<lower=0> a; // zero-point as fixed (e.g. 1)
}

transformed data {
  real eps = 1e-08; // small value for computational stability (avoid exact 0 and 1 in certain cases)
  real C = 1000; // maximum value for rho
}

parameters {
  // coefficients for beta distribution mean
  vector[n_var/2] beta_1; // first order coefficients
  vector<upper=0>[n_var/2] beta_2;

  // coefficients for probability of suitability
  vector[n_var/2] beta_pi_1;
  vector<upper=0>[n_var/2] beta_pi_2;

  // coefficients for modeling variance of Beta
  vector[n_var/2] beta_rho_1;
  vector[n_var/2] beta_rho_2;
  
  real alpha;
  real alpha_pi; // intercept for probability of suitability
  real alpha_rho; // intercept for rho
}

transformed parameters {
  // mean for beta distribution
  vector[N] mu;
  mu = inv_logit(alpha + X*append_row(beta_1,beta_2));
  
  // probability of suitability
  vector[N] prob_suit;
  prob_suit = inv_logit(alpha_pi + X*append_row(beta_pi_1,beta_pi_2));
  
  // rho parameters
  vector[N] rho;
  rho = C*inv_logit(alpha_rho + X*append_row(beta_rho_1,beta_rho_2));
}

model {
  // priors for coefficients
  beta_1 ~ normal(0,sqrt(2));
  beta_2 ~ normal(0,sqrt(2));
  beta_pi_1 ~ normal(0,sqrt(2));
  beta_pi_2 ~ normal(0,sqrt(2));
  beta_rho_1 ~ normal(0,sqrt(2));
  beta_rho_2 ~ normal(0,sqrt(2));
  
  alpha ~ normal(0,sqrt(2));
  alpha_pi ~ normal(-1,sqrt(0.25));
  alpha_rho ~ normal(0,sqrt(4));

  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter for beta distribution
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,eps), 1 - eps); // restrict to open (0,1)
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = prob_suit[i];
    pi_i = fmin(fmax(pi_i,eps), 1 - eps); // restrict to open (0,1)
    
    // calculate rho parameter
    real rho_i;
    rho_i = rho[i];
    rho_i = fmax(eps,rho_i);

    // likelihood terms
    if (y[i] == 0) {
      // logSumExp as in https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      //target += log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i));
      target += log_sum_exp(log(1-pi_i), log(pi_i) + log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i))));
    
    // since Beta-distribution has support on open interval (0,1), take y = 1 cases also separately, substract a small value
    } else if (y[i]==1) { 
      (y[i]-eps+a)/(a+1) ~ beta(mu_i*rho_i, (1-mu_i)*rho_i);
      target += log(pi_i);
      
    // y \in (0,1)
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho_i, (1-mu_i)*rho_i);
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
    
    real rho_i;
    rho_i = rho[i];
    rho_i = fmax(eps,rho_i);
    //rho_i = fmax(1e-08,exp(fmin(logit_rho[i],8)));
    
    // likelihood terms
    if (y[i] == 0) {
      //log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i));
      log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i))));
    } else if (y[i] == 1) {
      log_lik[i] = beta_lpdf((y[i]-eps+a)/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i) - log(1+a) + log(pi_i);
    } else {
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho_i, (1-mu_i)*rho_i) - log(1+a) + log(pi_i); // scale to [0,1] interval
    }
  }
}


