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
      // y_scaled = (y-eps+a)/(1+a); // this is how I scaled in thesis
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
  vector<lower=0,upper=1>[N] y; // species percentages
  matrix[N,n_var] X; // covariate matrix
  matrix[n_obs_grid,2] s; // coordinates of observed grid cells
  real<lower=0> a; // zero-point as fixed (e.g. 1)
  matrix[N,n_obs_grid] P; // binary matrix (nxm) telling to which grid cell each data point belongs 
}

transformed data {
  real eps = 1e-08; // small value for computational stability (avoid exact 0 and 1 in certain cases)
  array[n_obs_grid] vector[2] s_array; // kernel function takes input as array of vectors
  for(i in 1:n_obs_grid) {
    s_array[i] = to_vector(s[i]);
  }
}

parameters {
  // coefficients for beta distribution mean
  vector[n_var/2] beta_1;
  vector<upper=0>[n_var/2] beta_2;
  
  real alpha;
  real<lower=0> rho; // scale parameter for beta distribution
  
  real<lower=0> l; // length-scale parameter
  real<lower=0> s2_cf; // magnitude of the covariance function
  vector[n_obs_grid] z; // N(0,1) for creating the random effects
  //vector[n_obs_grid] phi; // spatial effects modeled excplicitely
}

transformed parameters {
  vector[n_obs_grid] phi;// random effects
  vector[N] mu;

  {
    matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
    matrix[n_obs_grid,n_obs_grid] L; // Cholesky of covariance matrix

    // construct K
    K = gp_exponential_cov(s_array,sqrt(s2_cf),l);

    // for (i in 1:n_obs_grid) {
    //   for (j in i:n_obs_grid) {
    //     K[i,j] = s2_cf*exp(-(distance(s[i],s[j]))/l);
    //     K[j,i] = s2_cf*exp(-(distance(s[i],s[j]))/l);
    //   }
    // }

    K = K + diag_matrix(rep_vector(eps,n_obs_grid)); // jitter for stability
    
    // generate the random effects
    L = cholesky_decompose(K);
    phi = L*z; // follows N(0,K)
    
    // mean parameter for the beta
    mu = inv_logit(alpha + X * append_row(beta_1,beta_2) + P*phi);
  }
}

model {
  // spatial random effect
  // matrix[n_obs_grid,n_obs_grid] K; // covariance matrix
  // //matrix[n_obs_grid,n_obs_grid] L; // Cholesky of covariance matrix
  // 
  // // construct K that depends on s2, l
  // for (i in 1:n_obs_grid) {
  //   for (j in i:n_obs_grid) {
  //     K[i,j] = s2_cf*exp(-(distance(s[i],s[j]))/l);
  //     K[j,i] = s2_cf*exp(-(distance(s[i],s[j]))/l);
  //   }
  // }
  
  //K = K + diag_matrix(rep_vector(1e-08,n_obs_grid)); // add jitter for stability
  //L = cholesky_decompose(K);
  //phi ~ multi_normal_cholesky(rep_vector(0,n_obs_grid),L); //GP(0,K)
  //phi ~ multi_normal(rep_vector(0,n_obs_grid),K); //GP(0,K)
  
  // priors for coefficients
  beta_1 ~ normal(0,sqrt(1));
  beta_2 ~ normal(0,sqrt(1));

  alpha ~ normal(0,sqrt(2)); 
  rho ~ cauchy(0,sqrt(10));
  
  z ~ normal(0,1);
  
  // priors for covariance function
  s2_cf ~ student_t(4,0,sqrt(0.1));
  l ~ cauchy(8,sqrt(50));
  
  for (i in 1:N) {
    target += lc_beta_density(y[i],mu[i],rho,a,eps,N);
    // // calculate the mean parameter for beta distribution
    // real mu_i;
    // mu_i = mu[i];
    // mu_i = fmin(fmax(mu_i,eps), 1 - eps); // restrict to open (0,1)
    // 
    // if (y[i] == 0) {
    //   // target += beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho);
    //   target += log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho)));
    //   
    // // since Beta-distribution has support on open interval (0,1), take y = 1 cases also separately, substract a small value
    // } else if (y[i] ==  1) {
    //   // substract a small value
    //   (y[i]-eps+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho);
    // 
    // // y \in (0,1)  
    // } else {
    //   (y[i]+a)/(a+1) ~ beta(mu_i*rho, (1-mu_i)*rho); // scale to [0,1] interval
    // }
  }
}

generated quantities {
  // log-likelihood for LOO-calculations
  vector[N] log_lik;
  
  for (i in 1:N) {
    log_lik[i] = lc_beta_density(y[i],mu[i],rho,a,eps,N);
    // real mu_i;
    // mu_i = mu[i];
    // mu_i = fmin(fmax(mu_i,1e-08), 1 - eps);
    // 
    // if (y[i] == 0) {
    //   //log_lik[i] = beta_lcdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho);
    //   log_lik[i] = log(fmax(1e-20, beta_cdf(a/(a+1) | mu_i*rho, (1-mu_i)*rho)));
    // } else if (y[i]==1) {
    //   log_lik[i] = beta_lpdf((y[i]-eps+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a);
    // } else {
    //   log_lik[i] = beta_lpdf((y[i]+a)/(a+1) | mu_i*rho, (1-mu_i)*rho) - log(1+a); // scale to [0,1] interval
    // }
  }
}
