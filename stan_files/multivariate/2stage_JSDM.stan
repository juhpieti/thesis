
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of data points
  int<lower=0> n_var; // number of variables
  int<lower=0> J; // number of species
  matrix<lower=0,upper=100>[N,J] Y; // species percentages
  vector<lower=0>[N] y_sum; // total percent cover (sum(y))
  matrix[N,n_var] X; // covariate matrix
}

transformed data {
  real eps = 1e-12; // for y = 0
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[n_var/2] beta_1; // first order terms
  vector<upper=0>[n_var/2] beta_2;
  real alpha;
  
  real<lower=0> s; // standard deviation
}

transformed parameters{
  vector[N] mu;
  mu = alpha + X*append_row(beta_1,beta_2);
}

model {
  log(y_sum + eps) ~ normal(mu, s);
  
  beta_1 ~ normal(0,sqrt(2));
  beta_2 ~ normal(0,sqrt(2));
  
  alpha ~ normal(0,10);
  
  //s ~ cauchy(0,2);
  s~student_t(4,0,1);
}

