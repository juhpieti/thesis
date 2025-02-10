### two different species, three models, differences?
### - 1) in loo
### - 2) post predictive checks
### - 3) coefficients...? 

# load in packages
library(terra)
library(loo)
library(ggplot2)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in helper functions
source("codes/helpers.R")

# load in the training data
#load("data/estonia_sub/estonia_sub_df.RData")
load("data/estonia_new/train_2020_2021.Rdata")

estonia_sub <- df_sub

# in the data there are three species with less than 5 observations, I will delete them
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- estonia_sub[,11:19]
X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

# for saving leave-one-out CVs
loo_amphi <- c()

### 1) Amphibalanus

amphi <- estonia_sub[,"Amphibalanus improvisus"]
amphi.01 <- amphi/100 # beta regression takes values from [0,1]

### RUN THE MODELS

# # 1) binomial regression
# amphi_dat.bin <- list(N = nrow(X.sec_ord),
#                       n_var = ncol(X.sec_ord),
#                       y = amphi,
#                       X = X.sec_ord)
# 
# mod_amphi.bin <- stan("stan_files/binomial_regression_iid_noise.stan", data = amphi_dat.bin, chains = 4, iter = 500, seed = 42)
# 
# # save loo
# loo_amphi <- c(loo_amphi, loo(mod_amphi.bin)$elpd_loo)
# 
# # posterior predictive
# par(mfrow = c(4,4),
#     mar = c(2,4,2,0))
# hist(amphi, breaks = 10, xlim = c(0,100), main = "obs", ylim = c(0,500))
# 
# alpha.sam <- as.matrix(mod_amphi.bin, pars = c("alpha"))
# beta.sam <- as.matrix(mod_amphi.bin, pars = c("beta"))
# s2.sam <- as.matrix(mod_amphi.bin, pars = c("s2"))
# 
# n_rep <- 15
# rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
# for (i in 1:n_rep) {
#   idx <- rep_idx[i]
#   y_rep <- c()
#   for (j in 1:length(amphi)) {
#     eps <- rnorm(1,0,sqrt(s2.sam[idx,]))
#     #eps <- 0
#     mu_i <- inv_logit(alpha.sam[idx, ] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,] + eps)
#     y_rep <- c(y_rep,rbinom(1,100,mu_i))
#   }
#   hist(y_rep, breaks = 10, xlim = c(0,100), main = paste0("rep",i), ylim = c(0,500))
# }
# 
# ### visualize
# plot_responses_binomial_regression(mod_amphi.bin,X.sec_ord,X,TRUE,TRUE,0,60,0,60,3,4)
# 
# ### in average conditions, distribution of y?
# average_distribution_binomial_regression(mod_amphi.bin,"amphi")

# 2) BETA REGRESSION
amphi_dat.beta <- list(N = nrow(X.sec_ord),
                 n_var = ncol(X.sec_ord),
                 y = amphi.01,
                 X = X.sec_ord,
                 a = 1)

mod_amphi.beta <- stan("stan_files/left_censored_beta_regression.stan", data = amphi_dat.beta, chains = 4, iter = 250, seed = 42)

# save loo value
loo_amphi <- c(loo_amphi,loo(mod_amphi.beta)$elpd_loo)

### posterior predictive check
# observed coverages
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(amphi.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))


alpha.sam <- as.matrix(mod_amphi.beta, pars = c("alpha"))
beta.sam <- as.matrix(mod_amphi.beta, pars = c("beta"))
rho.sam <- as.matrix(mod_amphi.beta, pars = c("rho"))

a <- 1

# make n_rep replications of new data, compare to observed coverages
n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE) #randomly take 15 sets of parameters
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
    rho_i <- rho.sam[idx,]
    V_i <- rbeta(1,mu_i*rho_i,(1-mu_i)*rho_i)
    y_rep_i <- max(0,(a+1)*V_i - a)
    y_rep <- c(y_rep, y_rep_i)
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_beta_regression(mod_amphi.beta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"expectation")
plot_responses_beta_regression(mod_amphi.beta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"probzero")

### plot distribution of y in average conditions (X = 0)
average_distribution_beta_regression(mod_amphi.beta,1,"amphi")

### try new predict function

# 3) ZI-Beta Regression

mod_amphi_ZIbeta <- stan("stan_files/zero_inflated_left_censored_beta_regression.stan",data = amphi_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_amphi <- c(loo_amphi,loo(mod_amphi_ZIbeta)$elpd_loo)

### posterior predictive checks
# histogram of coverages in observed data
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(amphi.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("alpha"))
beta.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("beta"))
rho.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("rho"))
alpha_pi.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("alpha_pi"))
beta_pi.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("beta_pi"))

a <- 1

# make n_rep replications of new data, compare to observed coverages
n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    # probability of suitability
    pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,])
    Z_i <- rbinom(1,1,pi_i) 
    if (Z_i == 0) { # location non-suitable => y=0
      y_rep <- c(y_rep,0)
    } else { # location suitable, sample value
      mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
      rho_i <- rho.sam[idx,]
      V_i <- rbeta(1,mu_i*rho_i,(1-mu_i)*rho_i)
      y_rep_i <- max(0,(a+1)*V_i - a)
      y_rep <- c(y_rep, y_rep_i)
    }
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### THERE WAS SOME PROBLEM SUCH THAT PARAMETERS WERE STUCK AT ZERO VALUES?
# y_rep <- c()
# for (j in 1:nrow(X.sec_ord)) {
#   pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,])
#   Z_i <- rbinom(1,1,pi_i)
#   if (Z_i == 0) {
#     y_rep <- c(y_rep,0)
#   } else {
#     mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
#     rho_i <- rho.sam[idx,]
#     V_i <- rbeta(1,mu_i*rho_i,(1-mu_i)*rho_i)
#     y_rep_i <- max(0,(a+1)*V_i - a)
#     y_rep <- c(y_rep, y_rep_i)
#   }
# }

# plot different sort of responses, expectation and probability of zero
plot_responses_ZI_beta_regression(mod_amphi_ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,type="expectation")
plot_responses_ZI_beta_regression(mod_amphi_ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,type="probzero")

average_distribution_ZI_beta_regression(mod_amphi_ZIbeta,1,"amphi")

### responses separately for suitability and abundance
plot_separate_responses_ZI_beta_regression(mod_amphi_ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1)

# save leave-one-out CV values
save(loo_amphi, file = "results/demo/amphi/loo_amphi.Rdata")


### Chara aspera

loo_chara <- c()

chara <- estonia_sub[,"Chara aspera"]
chara.01 <- chara/100 #beta regression takes values in [0,1]

### run models

# # 1) binomial regression
# chara_dat.bin <- list(N = nrow(X.sec_ord),
#                      n_var = ncol(X.sec_ord),
#                      y = chara,
#                      X = X.sec_ord)
# 
# mod_chara.bin <- stan("stan_files/binomial_regression_iid_noise.stan", data = chara_dat.bin, chains = 4, iter = 500, seed = 42)
# 
# # save LOO values
# loo_chara <- c(loo_chara, loo(mod_chara.bin)$elpd_loo)
# 
# # posterior predictive
# par(mfrow = c(4,4),
#     mar = c(2,4,2,0))
# hist(chara, breaks = 10, xlim = c(0,100), main = "obs", ylim = c(0,500))
# 
# alpha.sam <- as.matrix(mod_chara.bin, pars = c("alpha"))
# beta.sam <- as.matrix(mod_chara.bin, pars = c("beta"))
# s2.sam <- as.matrix(mod_chara.bin, pars = c("s2"))
# 
# n_rep <- 15
# rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
# for (i in 1:n_rep) {
#   idx <- rep_idx[i]
#   y_rep <- c()
#   for (j in 1:length(amphi)) {
#     eps <- rnorm(1,0,sqrt(s2.sam[idx,]))
#     #eps <- 0
#     mu_i <- inv_logit(alpha.sam[idx, ] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,] + eps)
#     y_rep <- c(y_rep,rbinom(1,100,mu_i))
#   }
#   hist(y_rep, breaks = 10, xlim = c(0,100), main = paste0("rep",i), ylim = c(0,500))
# }
# 
# ### visualize
# plot_responses_binomial_regression(mod_chara.bin,X.sec_ord,X,TRUE,TRUE,0,60,0,50,3,4)
# 
# ### in average conditions, distribution of y?
# average_distribution_binomial_regression(mod_chara.bin,"amphi")

# 2) Beta Regression
chara_dat.beta <- list(N = nrow(X.sec_ord),
                       n_var = ncol(X.sec_ord),
                       y = chara.01,
                       X = X.sec_ord,
                       a = 1)

mod_chara.beta <- stan("stan_files/left_censored_beta_regression.stan", data = chara_dat.beta, chains = 4, iter = 250, seed = 42)

# save LOO values
loo_chara <- c(loo_chara,loo(mod_chara.beta)$elpd_loo)

### posterior predictive checks
# histogram of coverages from real data
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(chara.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_chara.beta, pars = c("alpha"))
beta.sam <- as.matrix(mod_chara.beta, pars = c("beta"))
phi.sam <- as.matrix(mod_chara.beta, pars = c("phi"))

a <- 1

# make n_rep replications of new data, compare to observed coverages
n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
    phi_i <- phi.sam[idx,]
    V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
    y_rep_i <- max(0,(a+1)*V_i - a)
    y_rep <- c(y_rep, y_rep_i)
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_beta_regression(mod_chara.beta,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1)

### plot distribution of y in average conditions (X=0)
average_distribution_beta_regression(mod_chara.beta,1,"chara")

# 3) ZI-Beta Regression

mod_chara.ZIbeta <- stan("stan_files/zero_inflated_left_censored_beta_regression.stan",data = chara_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_chara <- c(loo_chara,loo(mod_chara.ZIbeta)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(chara.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_chara.ZIbeta, pars = c("alpha"))
beta.sam <- as.matrix(mod_chara.ZIbeta, pars = c("beta"))
phi.sam <- as.matrix(mod_chara.ZIbeta, pars = c("phi"))
alpha_pi.sam <- as.matrix(mod_chara.ZIbeta, pars = c("alpha_pi"))
beta_pi.sam <- as.matrix(mod_chara.ZIbeta, pars = c("beta_pi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,])
    Z_i <- rbinom(1,1,pi_i)
    if (Z_i == 0) {
      y_rep <- c(y_rep,0)
    } else {
      mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
      phi_i <- phi.sam[idx,]
      V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
      y_rep_i <- max(0,(a+1)*V_i - a)
      y_rep <- c(y_rep, y_rep_i)
    }
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

plot_responses_ZI_beta_regression(mod_chara.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)
average_distribution_ZI_beta_regression(mod_chara.ZIbeta,1,"chara")

### responses separately for suitability and abundance
plot_separate_responses_ZI_beta_regression(mod_chara.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)

save(loo_chara, file = "results/demo/chara/loo_chara.Rdata")


### save the coeffs
# coeffs_amphi <- c()
# coeffs_chara <- c()
# 
# coeffs_amphi <- cbind(get_posterior_mean(mod_amphi.bin, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_amphi.beta, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_amphi_ZIbeta, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_amphi_ZIbeta, pars = c("alpha_pi","beta_pi"))[,5])
# 
# colnames(coeffs_amphi) <- c("Bin","Beta","ZIBeta","ZIBeta (pi)")
# rownames(coeffs_amphi) <- c("intcpt",colnames(X.sec_ord))
# save(coeffs_amphi, file = "results/demo/amphi/coeffs_amphi.Rdata")
# 
# 
# coeffs_chara <- cbind(get_posterior_mean(mod_chara.bin, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_chara.beta, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_chara.ZIbeta, pars = c("alpha","beta"))[,5],
#                       get_posterior_mean(mod_chara.ZIbeta, pars = c("alpha_pi","beta_pi"))[,5])
# 
# colnames(coeffs_chara) <- c("Bin","Beta","ZIBeta","ZIBeta (pi)")
# rownames(coeffs_chara) <- c("intcpt",colnames(X.sec_ord))
# save(coeffs_amphi, file = "results/demo/amphi/coeffs_amphi.Rdata")


### 3) Mytilus trossulus

loo_mytilus <- c()

mytilus <- estonia_sub[,"Mytilus trossulus"]
mytilus.01 <- mytilus/100

# 2) beta regression
mytilus_dat.beta <- list(N = nrow(X.sec_ord),
                       n_var = ncol(X.sec_ord),
                       y = mytilus.01,
                       X = X.sec_ord,
                       a = 1)

mod_mytilus.beta <- stan("stan_files/left_censored_beta_regression.stan", data = mytilus_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_mytilus <- c(loo_mytilus,loo(mod_mytilus.beta)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(mytilus.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_mytilus.beta, pars = c("alpha"))
beta.sam <- as.matrix(mod_mytilus.beta, pars = c("beta"))
phi.sam <- as.matrix(mod_mytilus.beta, pars = c("phi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
    phi_i <- phi.sam[idx,]
    V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
    y_rep_i <- max(0,(a+1)*V_i - a)
    y_rep <- c(y_rep, y_rep_i)
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_beta_regression(mod_mytilus.beta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)

### plot distribution of y in average conditions
average_distribution_beta_regression(mod_mytilus.beta,1,"chara")

# 3) zero-inflated beta regression

mod_mytilus.ZIbeta <- stan("stan_files/zero_inflated_left_censored_beta_regression.stan",data = mytilus_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_mytilus <- c(loo_mytilus,loo(mod_mytilus.ZIbeta)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(mytilus.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_mytilus.ZIbeta, pars = c("alpha"))
beta.sam <- as.matrix(mod_mytilus.ZIbeta, pars = c("beta"))
phi.sam <- as.matrix(mod_mytilus.ZIbeta, pars = c("phi"))
alpha_pi.sam <- as.matrix(mod_mytilus.ZIbeta, pars = c("alpha_pi"))
beta_pi.sam <- as.matrix(mod_mytilus.ZIbeta, pars = c("beta_pi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,])
    Z_i <- rbinom(1,1,pi_i)
    if (Z_i == 0) {
      y_rep <- c(y_rep,0)
    } else {
      mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
      phi_i <- phi.sam[idx,]
      V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
      y_rep_i <- max(0,(a+1)*V_i - a)
      y_rep <- c(y_rep, y_rep_i)
    }
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

plot_responses_ZI_beta_regression(mod_mytilus.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)
average_distribution_ZI_beta_regression(mod_mytilus.ZIbeta,1,"chara")

### responses separately for suitability and abundance
plot_separate_responses_ZI_beta_regression(mod_mytilus.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)

save(loo_mytilus, file = "results/demo/mytilus/loo_mytilus.Rdata")


### 4) Zostera Marina

loo_zostera <- c()

zostera <- estonia_sub[,"Zostera marina"]
zostera.01 <- zostera/100

# 2) beta regression
zostera_dat.beta <- list(N = nrow(X.sec_ord),
                         n_var = ncol(X.sec_ord),
                         y = zostera.01,
                         X = X.sec_ord,
                         a = 1)

mod_zostera.beta <- stan("stan_files/left_censored_beta_regression.stan", data = zostera_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_zostera <- c(loo_zostera,loo(mod_zostera.beta)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(zostera.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_zostera.beta, pars = c("alpha"))
beta.sam <- as.matrix(mod_zostera.beta, pars = c("beta"))
phi.sam <- as.matrix(mod_zostera.beta, pars = c("phi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
    phi_i <- phi.sam[idx,]
    V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
    y_rep_i <- max(0,(a+1)*V_i - a)
    y_rep <- c(y_rep, y_rep_i)
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_beta_regression(mod_zostera.beta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)

### plot distribution of y in average conditions
average_distribution_beta_regression(mod_zostera.beta,1,"chara")

# 3) zero-inflated beta regression

mod_zostera.ZIbeta <- stan("stan_files/zero_inflated_left_censored_beta_regression.stan",data = zostera_dat.beta, chains = 4, iter = 250, seed = 42)

#save loo
loo_zostera <- c(loo_zostera,loo(mod_zostera.ZIbeta)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(zostera.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_zostera.ZIbeta, pars = c("alpha"))
beta.sam <- as.matrix(mod_zostera.ZIbeta, pars = c("beta"))
phi.sam <- as.matrix(mod_zostera.ZIbeta, pars = c("phi"))
alpha_pi.sam <- as.matrix(mod_zostera.ZIbeta, pars = c("alpha_pi"))
beta_pi.sam <- as.matrix(mod_zostera.ZIbeta, pars = c("beta_pi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,])
    Z_i <- rbinom(1,1,pi_i)
    if (Z_i == 0) {
      y_rep <- c(y_rep,0)
    } else {
      mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,])
      phi_i <- phi.sam[idx,]
      V_i <- rbeta(1,mu_i*phi_i,(1-mu_i)*phi_i)
      y_rep_i <- max(0,(a+1)*V_i - a)
      y_rep <- c(y_rep, y_rep_i)
    }
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

plot_responses_ZI_beta_regression(mod_zostera.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)
average_distribution_ZI_beta_regression(mod_zostera.ZIbeta,1,"chara")

### responses separately for suitability and abundance
plot_separate_responses_ZI_beta_regression(mod_zostera.ZIbeta,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)

save(loo_zostera, file = "results/demo/zostera/loo_zostera.Rdata")


### save the loo table
loo_table <- rbind(loo_amphi, loo_chara, loo_mytilus, loo_zostera)
rownames(loo_table) <- c("Amphibalanus improvisus","Chara aspera","Mytilus trossulus","Zostera marina")
colnames(loo_table) <- c("Beta","ZI-Beta")
save(loo_table, file = "results/demo/loo_table.Rdata")


####################################################################################################################################
### test area

# clayplate = 20, others are in the mean
# x_pred <- colMeans(X)
# x_pred["clayplate"] <- 20
# x_pred <- matrix(x_pred, ncol = length(x_pred))
# x_pred <- as.data.frame(x_pred) 
# colnames(x_pred) <- colnames(X)
# x_pred <- scale_covariates(X,x_pred)
# x_pred <- add_second_order_terms(x_pred, colnames(X))
# 
# 
# ### mod 1
# 
# alpha.sam <- as.matrix(mod_amphi.bin, pars = c("alpha"))
# beta.sam <- as.matrix(mod_amphi.bin, pars = c("beta"))
# s2.sam <- as.matrix(mod_amphi.bin, pars = c("s2"))
# 
# alpha_mean <- mean(alpha.sam)
# beta_mean <- colMeans(beta.sam)
# s2_mean <- mean(s2.sam)
# 
# n_rep <- 100
# y_rep <- rbinom(n_rep,100,inv_logit(alpha_mean + as.numeric(x_pred[1,]) %*% beta_mean + rnorm(1,0,sqrt(s2_mean))))
# par(mfrow = c(1,1))
# hist(y_rep)
# mean(y_rep)
# 
# ### MOD 2
# 
# alpha.sam <- as.matrix(mod_amphi.beta, pars = c("alpha"))
# beta.sam <- as.matrix(mod_amphi.beta, pars = c("beta"))
# phi.sam <- as.matrix(mod_amphi.beta, pars = c("phi"))
# 
# alpha_mean <- mean(alpha.sam)
# beta_mean <- colMeans(beta.sam)
# phi_mean <- mean(phi)
# 
# mu <- inv_logit(alpha_mean + as.numeric(x_pred[1,]) %*% beta_mean )
# 
# V_rep <- rbeta(100,mu*phi_mean,(1-mu)*phi_mean)
# y_rep <- c()
# for (i in 1:length(V_rep)) {
#   y_rep <- c(y_rep, max(0,V_rep[i]*(a+1)-a))
# }
# 
# hist(y_rep)
# mean(y_rep)
# 
# 
# grid_01 <- seq(0.01,1-0.01,0.01)
# 
# sum(0.01*dbeta((grid_01+a)/(a+1),mu*phi_mean,(1-mu)*phi_mean)/(a+1))
# 
# 
# ### mod 3
# 
# alpha.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("alpha"))
# alpha_pi.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("alpha_pi"))
# beta.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("beta"))
# beta_pi.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("beta_pi"))
# phi.sam <- as.matrix(mod_amphi_ZIbeta, pars = c("phi"))
# 
# alpha_mean <- mean(alpha.sam)
# alpha_pi_mean <- mean(alpha_pi.sam)
# beta_mean <- colMeans(beta.sam)
# beta_pi_mean <- colMeans(beta_pi.sam)
# phi_mean <- mean(phi)
# 
# mu <- inv_logit(alpha_mean -0.8827*beta_mean[1])
# prob_suit <- inv_logit(alpha_pi_mean - 0.8827*beta_pi_mean[1])
# 
# calc_dens_test <- function(x) (x*prob_suit*dbeta((x+a)/(a+1),mu*phi_mean,(1-mu)*phi_mean)/(a+1))
# integrate(calc_dens_test,0,1)
