//
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  real<lower=0> a; // zero-point as fixed (e.g. 1)
}

parameters {
  vector[n_var/2] beta_1; // first order terms
  //vector<upper=0>[n_var/2] beta_2;
  vector[n_var/2] logneg_beta_2; //unconstrained version, theta = log(-B) <=> B = -exp(theta)
  real alpha;
  real<lower=0> rho; // scale parameter for beta distribution
}

transformed parameters{
  vector[N] mu;
  mu = inv_logit(alpha + X*append_row(beta_1,-exp(logneg_beta_2)));
}

model {
  // priors for coefficients
  //beta_1 ~ multi_normal(rep_vector(0,n_var/2),diag_matrix(rep_vector(sqrt(2),n_var/2)));
  beta_1 ~ normal(0,sqrt(2));
  //beta_2 ~ multi_normal(rep_vector(0,n_var/2),diag_matrix(rep_vector(sqrt(2),n_var/2)));
  logneg_beta_2 ~ normal(0,1); 

  alpha ~ normal(0,sqrt(2)); // ORIGINAL
  // rho ~ inv_gamma(0.5,0.5);
  rho ~ cauchy(0,sqrt(5));
  
  // likelihood
  for (i in 1:N) {
    // calculate the mean parameter
    real mu_i;
    //mu_i = inv_logit(alpha + X[i,1:(n_var/2)] * beta_1 + X[i,(1+n_var/2):n_var] * beta_2);
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
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
  // return beta_2 in constrained space
  vector[n_var/2] beta_2;
  beta_2 = -exp(logneg_beta_2);
  
  vector[N] log_lik;

  for (i in 1:N) {
    real mu_i;
    //mu_i = inv_logit(alpha + X[i,1:(n_var/2)] * beta_1 + X[i,(1+n_var/2):n_var] * beta_2);
    mu_i = mu[i];
    mu_i = fmin(fmax(mu_i,1e-08), 1 - 1e-08);
    
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
