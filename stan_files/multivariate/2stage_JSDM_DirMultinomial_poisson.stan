
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  int<lower=0> J; // number of species
  int<lower=0> n_f; // number of latent factors
  //matrix<lower=0,upper=100>[N,J] Y; // species percentages
  array[N, J] int<lower=0, upper=100> Y; // species percentages
  //vector<lower=0>[N] y_sum; // total percent cover (sum(y))
  array[N] int<lower=0> y_sum; // total percent cover (sum(y))
  matrix[N,n_var] X; // covariate matrix
}

transformed data {
  real eps = 1e-12; // for y = 0
}

parameters {
  // coefficients for Dirichlet distribution
  matrix[n_var/2,J] beta_1; // first order terms
  matrix<upper=0>[n_var/2,J] beta_2; // second order terms, restricted negative (bell-shaped curves based on ecological niche theory)
  vector[J] alpha; // intercept terms
  real<lower=0> scale_Dir; // scaling factor
  
  // latent factors & species loadings
  matrix[N,n_f] Z; // n_f latent factors for each sampling location
  matrix[n_f,J] Lambda; // species loadings
  
  // parameters for the total coverage model
  vector[n_var/2] beta_M_1; // first order terms
  vector<upper=0>[n_var/2] beta_M_2; // second order terms, restricted negative (bell-shaped curves based on ecological niche theory)
  real alpha_M; // intercept
  
  // Expected proportions
  // array[N] simplex[J] Pi; // N rows, each are J-vectors that sum up to 1 (simplex)
}

transformed parameters{
  // Parameters for the Dirichlet distribution
  matrix[N,J] Mu_Dir; // gather the mean parameters for beta distribution
  //array[N,J] real<lower=0,upper=1> Mu_Dir; // gather the mean parameters for beta distribution
  for (n in 1:N) {
    // Mu_Dir[n] = softmax(to_row_vector(alpha) + X[n]*append_row(beta_1,beta_2) + Z[n]*Lambda); // intercept + fixed effects + random effects
    Mu_Dir[n] = to_row_vector(softmax(to_vector(alpha) + to_vector(X[n]*append_row(beta_1,beta_2)) + to_vector(Z[n]*Lambda))); // intercept + fixed effects + random effects
  }
  
  // Parameters for the model for total coverage
  vector[N] mu_M;
  mu_M = exp(alpha_M + X*append_row(beta_M_1,beta_M_2));
}

model {
  // LIKELIHOODS
  for (i in 1:N) {
    // OBSERVATION MODEL
    //Y[i] ~ multinomial(to_vector(Pi[i,]));
    Y[i] ~ dirichlet_multinomial(scale_Dir*to_vector(Mu_Dir[i]));
    // PROPORTIONS
    //Pi[i,] ~ dirichlet(scale_Dir*Mu_Dir[i,]);
  }
  
  // TOTAL COVER
  y_sum ~ poisson(mu_M);
  
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
  
  // factor loadings
  to_vector(Z) ~ normal(0,1);
  //to_vector(Lambda) ~ normal(0,1);
  to_vector(Lambda) ~ normal(0,0.5);
}
