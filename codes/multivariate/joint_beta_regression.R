# load in packages
library(terra)
library(loo)
library(ggplot2)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in utility/helper functions
source("codes/helpers.R")

### DATA PREPARATION

# load in the dataset (n=100)
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
dim(train_n100)
colSums(train_n100[,20:38] > 0) #number of presences

# examine the covariates for correlation
library(corrplot)
X <- train_n100[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate the light level at the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# draw correlation plot
cor.mat <- cor(X)
corrplot(cor.mat, type = "upper", order = "hclust", method = "number")

# take a subset of covariates
cov_names <- c("depth","o2_bottom","bottomT","so_bottom","light_bottom")
X <- X[ ,cov_names]

# correlation plot for species percent covers
Y <- train_n100[,20:38]
Y <- Y[,colSums(Y>0)>12] #take only species with > 10 appearances

cor.mat <- cor(Y)
corrplot(cor.mat, type = "upper", order = "hclust", method = "number")

# take interesting, somewhat correlated species
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y <- Y[ ,sp_names]


### MODELING PART
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms

Y.scaled <- Y/100

data.list <- list(N = nrow(Y.scaled),
                  n_var = ncol(X.sec_ord),
                  J = ncol(Y.scaled),
                  n_f = 2,
                  Y = Y.scaled,
                  X = X.sec_ord,
                  a = 1)

n_chains <- 4
n_iter <- 2000

fit.beta.JSDM <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate.stan",
                      data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                      pars = c("Mu"), include = FALSE)

fit.beta.JSDM.rho.modeled <- stan("stan_files/multivariate/scaled_sigmoid/left_censored_beta_regression_multivariate.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)


### with hierarchical priors
fit.beta.JSDM.hier.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_hierarchical_priors.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)

### with gamma shrinkage prior for Lambda
fit.beta.JSDM.gamma.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_gamma_shrinkage_priors.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)

### SAVE THE MODEL FIT
subfolder <- paste0("n_",nrow(Y.scaled),"/",ncol(X),"_covariates/")

saveRDS(fit.beta.JSDM, file = paste0("models/multivariate/",subfolder,"JSDM_test.RDS"))
saveRDS(fit.beta.JSDM.rho.modeled, file = "models/multivariate/JSDM_test_rho_modeled.RDS")
saveRDS(fit.beta.JSDM.hier.priors, file = "models/multivariate/JSDM_test_hierarchical_priors.RDS")
saveRDS(fit.beta.JSDM.gamma.priors, file = "models/multivariate/JSDM_test_gamma_priors.RDS")




### EXAMINE THE MODEL FIT
fit.beta.JSDM <- readRDS("models/multivariate/JSDM_test.RDS")
loo(fit.beta.JSDM)

### SITE LOADINGS (color with some missing variables?)
Z.sam <- as.matrix(fit.beta.JSDM, pars = "Z")
Z.post_means <- colMeans(Z.sam)
Z.post_means <- matrix(Z.post_means, nrow = nrow(Y), byrow = FALSE)

plot(Z.post_means[,1],Z.post_means[,2], pch = 19, xlab = "latent factor 1", ylab = "latent factor 2")
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# Color points wrt. covariates not used in modeling
df <- as.data.frame(cbind(Z.post_means,train_n100[,c("no3_bottom","po4_bottom","current_bottom","chl_bottom")]))
colnames(df)[1:2] <- c("lf1","lf2")

plots <- list()
for (var in colnames(df)[3:6]) {
  plots[[var]] <- ggplot(data=df,aes(x = lf1, y = lf2)) + 
    geom_point(aes(colour=.data[[var]]))
}
gridExtra::grid.arrange(grobs = plots, ncol = 2)


### SPECIES LOADINGS
L.sam <- as.matrix(fit.beta.JSDM, pars = "Lambda")
L.post_means <- colMeans(L.sam)
L.post_means <- matrix(L.post_means, ncol = ncol(Y), byrow = FALSE)

plot(L.post_means[1,],L.post_means[2,], pch = 19, xlab = "loading 1", ylab="loading 2")
text(L.post_means[1,],L.post_means[2,], colnames(Y), cex = 0.8, pos = 3)
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# one can also calculate the residual covariance matrix
cor.mat <- cov2cor(t(L.post_means) %*% L.post_means)
colnames(cor.mat) <- rownames(cor.mat) <- colnames(Y)
cor.mat
