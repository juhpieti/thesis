#######################################################################################
### Script that plays with the beta distribution shapes wrt. mu,rho parametrization ###
#######################################################################################

### the effect of sample size (rho) with constant mu
mu <- 0.5
sample_size <- c(0.1,1,2,5,10,50)
x_grid <- seq(0,1,length=100)

cols <- rainbow(length(sample_size))
plot(NULL, xlim = 0:1, ylim = c(0,6))
for (i in 1:length(sample_size)) {
  ss <- sample_size[i]
  print(ss)
  lines(x_grid,dbeta(x_grid,mu*ss,(1-mu)*ss), type = "l", col = cols[i])
}
legend("topright", legend = paste0("sample size = ", sample_size), lty = 1, col = cols)

### Model with varying rho is implemented rho(x) = C*inverse_logit(a+xb)
### This parametrization requires fixing C (the maximum value for rho)
### Observed precent covers are reported as 1,5,10,15,...,95,100
### Considering the measurement error (5%), think about how large the rho parameter can be:

rho <- c(2,4,10,20,50,100,250,500,1000)
mu <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
x_seq <- seq(0,1,length=200)
for (mu_i in mu) {
  par(mfrow = c(3,3))
  for (rho_i in rho) {
    prob_within_5 <- pbeta(mu_i+0.025, mu_i*rho_i, (1-mu_i)*rho_i) - pbeta(mu_i-0.025, mu_i*rho_i, (1-mu_i)*rho_i)
    plot(x_seq,dbeta(x_seq,mu_i*rho_i,(1-mu_i)*rho_i), type = "l", main = paste0("mu=",mu_i,"rho=",rho_i))
    legend("topright",legend = prob_within_5)
    abline(v=0.475,col="red",lty=2)
    abline(v=0.525,col="red",lty=2)
    
  }
}

### Model with rho(x) as function of covariates resulted in weirdly U-shaped response curves
### examine why this happens
mod.beta <- readRDS("models/n_500/M1/Cladophora_glomerata.rds")
mod.beta_rho <- readRDS("models/scaled_sigmoid/n_500/M1/Cladophora_glomerata.rds")

### draw curves of mu and rho wrt. salinity (6th covariate)
# mu related parameters
alpha.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("alpha")))
beta1.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("beta_1[6]")))
beta2.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("beta_2[6]")))

# mu related parameters for model with common rho
alpha_base.mean <- colMeans(as.matrix(mod.beta, pars = c("alpha")))
beta1_base.mean <- colMeans(as.matrix(mod.beta, pars = c("beta_1[6]")))
beta2_base.mean <- colMeans(as.matrix(mod.beta, pars = c("beta_2[6]")))

# rho related parameters
alpha_rho.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("alpha_rho")))
beta1_rho.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("beta_rho_1[6]")))
beta2_rho.mean <- colMeans(as.matrix(mod.beta_rho, pars = c("beta_rho_2[6]")))

# mean value for rho in model of common rho
rho_base.mean <- colMeans(as.matrix(mod.beta, pars = c("rho")))

# load in the training data
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
X <- train_n500[,11:19]

# make a grid for salinity
so_grid <- seq(min(X$so_bottom),max(X$so_bottom),length=200)
so_grid_scaled <- (so_grid - mean(X$so_bottom))/ sd(X$so_bottom)

# make a grid for mu
mu_grid <- inv_logit(alpha.mean + so_grid_scaled*beta1.mean + beta2.mean*so_grid_scaled^2)
mu_base_grid <- inv_logit(alpha_base.mean + so_grid_scaled*beta1_base.mean + beta2_base.mean*so_grid_scaled^2)

# make a grid for rho
rho_grid <- 1000*inv_logit(alpha_rho.mean + beta1_rho.mean*so_grid_scaled + beta2_rho.mean*so_grid_scaled^2)

# plot curves (mu)
plot(so_grid, mu_grid, type = "l", ylim = c(0,0.5))
lines(so_grid, mu_base_grid, col = "red")
legend("bottomleft", legend = c("common rho","varying rho"), lty = c(1,1), col = c("red","black"))

# plot curves (rho)
plot(so_grid, rho_grid, type = "l", ylim = c(0,100))
abline(h = rho_base.mean, lty = 2, col = "red")
legend("topright", legend = c("common rho","varying rho"), lty = c(2,1), col = c("red","black"))

### take three cases and examine the shape of the latent Beta distribution
### 1) salinity = 3.3
### 2) salinity = 5.5
### 3) salinity = 8.5

sal1 <- 3.3
sal2 <- 5.5
sal3 <- 8.5

sal1.scaled <- (sal1-mean(X$so_bottom))/sd(X$so_bottom)
sal2.scaled <- (sal2-mean(X$so_bottom))/sd(X$so_bottom)
sal3.scaled <- (sal3-mean(X$so_bottom))/sd(X$so_bottom)

seq.01 <- seq(0,1,length=200)

mu1 <- inv_logit(alpha.mean + beta1.mean*sal1.scaled + beta2.mean*sal1.scaled^2)
rho1 <- 1000*inv_logit(alpha_rho.mean + beta1_rho.mean*sal1.scaled + beta2_rho.mean*sal1.scaled^2)

plot(seq.01, dbeta(seq.01, mu1*rho1, (1-mu1)*rho1), type = "l")
abline(v=0.5,lty=2,col="red")

### calculate expectation first by sampling
a <- 1
V.sam <- rbeta(10000,mu1*rho1,(1-mu1)*rho1)
W.sam <- (a+1)*V.sam - a
Y.sam <- sapply(W.sam, FUN = function(w) (max(0,w)))
mean(Y.sam) #0.3037

### then with my function (just to test it works)
source("codes/helpers.R")
integrate_LC_beta(mu1,rho1,1,0.5,FALSE) #0.3066

### same thing with the second set of mu,rho values
mu2 <- inv_logit(alpha.mean + beta1.mean*sal2.scaled + beta2.mean*sal2.scaled^2)
rho2 <- 1000*inv_logit(alpha_rho.mean + beta1_rho.mean*sal2.scaled + beta2_rho.mean*sal2.scaled^2)

# plot latent beta distribution
plot(seq.01, dbeta(seq.01, mu2*rho2, (1-mu2)*rho2), type = "l")
abline(v=0.5,lty=2,col="red")

### calculate expectations
a <- 1
V.sam <- rbeta(10000,mu2*rho2,(1-mu2)*rho2)
W.sam <- (a+1)*V.sam - a
Y.sam <- sapply(W.sam, FUN = function(w) (max(0,w)))
mean(Y.sam) #0.0078
integrate_LC_beta(mu2,rho2,1,0.5,FALSE) #0.0075

### and with the third set of mu, rho values
mu3 <- inv_logit(alpha.mean + beta1.mean*sal3.scaled + beta2.mean*sal3.scaled^2)
rho3 <- 1000*inv_logit(alpha_rho.mean + beta1_rho.mean*sal3.scaled + beta2_rho.mean*sal3.scaled^2)

# plot latent beta distribution
plot(seq.01, dbeta(seq.01, mu3*rho3, (1-mu3)*rho3), type = "l")
abline(v=0.5,lty=2,col="red")

# calculate expectations
a <- 1
V.sam <- rbeta(10000,mu3*rho3,(1-mu3)*rho3)
W.sam <- (a+1)*V.sam - a
Y.sam <- sapply(W.sam, FUN = function(w) (max(0,w)))
mean(Y.sam) #0.328
integrate_LC_beta(mu3,rho3,1,0.5,FALSE) #0.307
