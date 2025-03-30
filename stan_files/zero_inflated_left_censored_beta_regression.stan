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
  //vector[n_var/2] logneg_beta_2; // second order coefficients
  
  // coefficients for probability of suitability
  vector[n_var/2] beta_pi_1;
  vector<upper=0>[n_var/2] beta_pi_2;
  //vector[n_var/2] logneg_beta_pi_2; 
  
  real alpha;
  real alpha_pi; // intercept for probability of suitability
  real<lower=0> rho; // scale parameter for beta distribution
}

transformed parameters {
  // mean for beta distribution
  vector[N] mu;
  //mu = inv_logit(alpha + X*append_row(beta_1,-exp(logneg_beta_2)));
  mu = inv_logit(alpha + X*append_row(beta_1,beta_2));
  // probability of suitability
  vector[N] prob_suit;
  //prob_suit = inv_logit(alpha_pi + X*append_row(beta_pi_1,-exp(logneg_beta_pi_2)));
  prob_suit = inv_logit(alpha_pi + X*append_row(beta_pi_1,beta_pi_2));
}

model {
  // priors for coefficients
  beta_1 ~ normal(0,sqrt(2));
  beta_2 ~ normal(0,sqrt(2));
  //beta_1 ~ double_exponential(0,1); //Laplace
  //beta_2 ~ double_exponential(0,1);
  //logneg_beta_2 ~ normal(0,1);
  beta_pi_1 ~ normal(0,sqrt(2));
  beta_pi_2 ~ normal(0,sqrt(2));
  //beta_pi_1 ~ double_exponential(0,1);
  //beta_pi_2 ~ double_exponential(0,1);
  //logneg_beta_pi_2 ~ normal(0,1);
  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  alpha_pi ~ normal(0,sqrt(2)); //
  rho ~ cauchy(0,sqrt(10));
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter for beta distribution
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08); // there was problems of getting exact 0/1 with beta distribution calls
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = prob_suit[i];
    pi_i = fmin(fmax(pi_i,1e-08), 1 - 1e-08);
    
    // likelihood terms
    if (y[i] == 0) {
      // logSumExp as in https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
      target += log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
    } else if (y[i]==1) { // stan declares beta on open interval (0,1)
      (y[i]-1e-08+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho);
      target += log(pi_i);
    } else {
      (y[i]+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho); // scale to [0,1] interval
      target += log(pi_i);
    }
  }
}

generated quantities {
  // return beta_2 and beta_pi_2 in restricted spae
  //vector[n_var/2] beta_2;
  //beta_2 = -exp(logneg_beta_2);
  
  //vector[n_var/2] beta_pi_2;
  //beta_pi_2 = -exp(logneg_beta_pi_2);
  
  vector[N] log_lik;
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
    // calculate the probability of suitability
    real pi_i;
    pi_i = prob_suit[i];
    pi_i = fmin(fmax(pi_i,1e-08), 1 - 1e-08);
    
    // likelihood terms
    if (y[i] == 0) {
      log_lik[i] = log_sum_exp(log(1-pi_i), log(pi_i) + beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho));
    } else if (y[i] == 1) {
      log_lik[i] = beta_lpdf((y[i]-1e-08+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a) + log(pi_i);
      //log_lik[i] += log(pi_i);
    } else {
      // I suspect the log of jacobian should be added! log(1/a+1)=log(1/2)=-0.69...!
      log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a) + log(pi_i); // scale to [0,1] interval
      //log_lik[i] += log(pi_i);
    }
  }
}


