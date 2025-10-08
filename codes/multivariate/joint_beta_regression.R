################################################################
### SCRIPT TO RUN MULTIVARIATE LEFT-CENSORED BETA REGRESSION ###
################################################################

# load in packages
library(terra)
library(loo)
library(ggplot2)
library(corrplot)


# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in utility/helper functions
source("codes/helpers.R")

############################# DEMO WITH FEW (J=4) SPECIES ##############################

### DATA PREPARATION ###
# load in the dataset (n=100)
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
dim(train_n100)
colSums(train_n100[,20:38] > 0) #number of presences


# examine the covariates for correlation (possibly to drop highly correlated ones for analysis)
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

# correlation plot for species percent covers (to start with, take some correlating species for modeling task)
Y <- train_n100[,20:38]
Y <- Y[,colSums(Y>0)>12] #take only species with > 10 appearances

cor.mat <- cor(Y)
corrplot(cor.mat, type = "upper", order = "hclust", method = "number")

# choose interesting, somewhat correlated species
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y <- Y[ ,sp_names]

### MODELING PART ###
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms
Y.scaled <- Y/100 #scale y to (0,1) for beta regression

# prepare data for stan
data.list <- list(N = nrow(Y.scaled),
                  n_var = ncol(X.sec_ord),
                  J = ncol(Y.scaled),
                  n_f = 2,
                  Y = Y.scaled,
                  X = X.sec_ord,
                  a = 1)

# stan input parameters
n_chains <- 4
n_iter <- 2000

# fit the model
fit.beta.JSDM <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate.stan",
                      data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                      pars = c("Mu"), include = FALSE)

# fit the alternative model where rho(x) is modeled by covariates
fit.beta.JSDM.rho.modeled <- stan("stan_files/multivariate/scaled_sigmoid/left_censored_beta_regression_multivariate.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)


# fit alternative model with hierarchical priors for beta
fit.beta.JSDM.hier.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_hierarchical_priors.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)

# fit the alternative model with gamma shrinkage prior for Lambda (from Ovaskainen book)
fit.beta.JSDM.gamma.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_gamma_shrinkage_priors.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)

### SAVE THE MODEL FIT(s) ###
subfolder <- paste0("n_",nrow(Y.scaled),"/",ncol(X),"_covariates/")

saveRDS(fit.beta.JSDM, file = paste0("models/multivariate/",subfolder,"JSDM_test.RDS"))
saveRDS(fit.beta.JSDM.rho.modeled, file = "models/multivariate/JSDM_test_rho_modeled.RDS")
saveRDS(fit.beta.JSDM.hier.priors, file = paste0("models/multivariate/",subfolder,"JSDM_test_hierarchical_priors.RDS"))
saveRDS(fit.beta.JSDM.gamma.priors, file = paste0("models/multivariate/",subfolder,"/JSDM_test_gamma_priors_nf4.RDS"))

### EXAMINE THE MODEL FIT ###
fit.beta.JSDM <- mod.JSDM_4species

fit.beta.JSDM <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test.RDS")
fit.beta.JSDM <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test_gamma_priors_nf2.RDS")
fit.beta.JSDM <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test_gamma_priors_nf4.RDS")


### VISUALIZE SITE LOADINGS (color with some of the missing variables?)
Z.sam <- as.matrix(fit.beta.JSDM, pars = "Z")
Z.post_means <- colMeans(Z.sam)
# reorder as matrix (N x n_latent_factors)
Z.post_means <- matrix(Z.post_means, nrow = nrow(Y), byrow = FALSE)

plot(Z.post_means[,1],Z.post_means[,2], pch = 19, xlab = "latent factor 1", ylab = "latent factor 2")
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# Color points wrt. covariates not used in modeling (to see if there is any pattern)
df <- as.data.frame(cbind(Z.post_means,train_n100[,c("no3_bottom","po4_bottom","current_bottom","chl_bottom")]))
colnames(df)[1:2] <- c("lf1","lf2")

plots <- list()
for (var in colnames(df)[3:6]) {
  plots[[var]] <- ggplot(data=df,aes(x = lf1, y = lf2)) + 
    geom_point(aes(colour=.data[[var]]))
}
gridExtra::grid.arrange(grobs = plots, ncol = 2)

### VISUALIZE SPECIES LOADINGS ###
L.sam <- as.matrix(fit.beta.JSDM, pars = "Lambda")
L.post_means <- colMeans(L.sam)
# reorder as matrix (n_latent_factors x n_species)
L.post_means <- matrix(L.post_means, ncol = ncol(Y), byrow = FALSE)

plot(L.post_means[1,],L.post_means[2,], pch = 19, xlab = "loading 1", ylab="loading 2")
text(L.post_means[1,],L.post_means[2,], colnames(Y), cex = 0.8, pos = 3)
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# one can also calculate the residual covariance matrix
cor.mat <- cov2cor(t(L.post_means) %*% L.post_means)
colnames(cor.mat) <- rownames(cor.mat) <- colnames(Y)
cor.mat


################################## TRY THE LIMITS OF THE MODEL BY INTRODUCING MORE SPECIES ################################


# load in the dataset (n=100)
load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
train <- train_n100_all_species

# load in the dataset (n=250)
load("data/estonia_new/train/train_2020_2021_all_species_n250.Rdata")
train <- train_n250_all_species

load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
train <- train_n500_all_species

dim(train)
colnames(train)
colSums(train[,20:71] > 0) #number of presences

# covariates
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate the light level at the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# correlation plot for species percent covers (to start with, take some correlating species for modeling task)
Y <- train[,20:71]
Y <- Y[,colSums(Y>0)>0] #take only species with > 0 appearances
cor.mat <- cor(Y)
par(mfrow = c(1,1))
corrplot(cor.mat, type = "upper", order = "hclust", method = "number", number.cex = 0.6)

### MODELING PART ###
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms
Y.scaled <- Y/100 #scale y to (0,1) for beta regression

# prepare data for stan
data.list <- list(N = nrow(Y.scaled),
                  n_var = ncol(X.sec_ord),
                  J = ncol(Y.scaled),
                  n_f = 2,
                  Y = Y.scaled,
                  X = X.sec_ord,
                  a = 1)

# stan input parameters
n_chains <- 4
n_iter <- 2000

# fit the model
fit.beta.JSDM <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate.stan",
                      data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                      pars = c("Mu"), include = FALSE)

n_chains <- 4
n_iter <- 500
start.time <- Sys.time()
fit.beta.JSDM.n500 <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate.stan",
                      data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                      pars = c("Mu"), include = FALSE)
end.time <- Sys.time()
print(end.time - start.time)

subfolder <- paste0("n_",nrow(Y.scaled),"/")
saveRDS(fit.beta.JSDM.n500, file = paste0("models/multivariate/",subfolder,"M1/JSDM_21_species.RDS"))

### SAVE THE MODEL FIT(s) ###
subfolder <- paste0("n_",nrow(Y.scaled),"/")

saveRDS(fit.beta.JSDM, file = paste0("models/multivariate/",subfolder,"M1/JSDM_21species_Lambda_N0_05.RDS"))

# # fit the alternative model where rho(x) is modeled by covariates
# fit.beta.JSDM.rho.modeled <- stan("stan_files/multivariate/scaled_sigmoid/left_censored_beta_regression_multivariate.stan",
#                                   data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                                   pars = c("Mu"), include = FALSE)

# 
# # fit alternative model with hierarchical priors for beta
fit.beta.JSDM.hier.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_hierarchical_priors.stan",
                                  data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                  pars = c("Mu"), include = FALSE)
# 
saveRDS(fit.beta.JSDM.hier.priors, file = paste0("models/multivariate/",subfolder,"M1/JSDM_21species_hierarchical_priors.RDS"))


 
# fit the alternative model with gamma shrinkage prior for Lambda (from Ovaskainen book)
fit.beta.JSDM.gamma.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_gamma_shrinkage_priors.stan",
                                   data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                   pars = c("Mu"), include = FALSE)

subfolder <- paste0("n_",nrow(Y.scaled),"/")
saveRDS(fit.beta.JSDM.gamma.priors, file = paste0("models/multivariate/",subfolder,"M1/JSDM_21_species_gamma_priors.RDS"))


# fit the alternative model with both gamma shrinkage for Lambda and hierarchical prior for beta
fit.beta.JSDM.hier.gamma.priors <- stan("stan_files/multivariate/left_censored_beta_regression_multivariate_gamma_shrinkage_and_hierarchical_priors.stan",
                                        data = data.list, chains = n_chains, iter = n_iter, seed = 42,
                                        pars = c("Mu"), include = FALSE)

subfolder <- paste0("n_",nrow(Y.scaled),"/")
saveRDS(fit.beta.JSDM.hier.gamma.priors, file = paste0("models/multivariate/",subfolder,"M1/JSDM_21_species_hier_and_gamma_priors.RDS"))


### examine the interspecific correlations
fit.beta.JSDM_21_species <- readRDS("models/multivariate/n_500/M1/JSDM_21_species.RDS")

### VISUALIZE SPECIES LOADINGS ###
L.sam <- as.matrix(fit.beta.JSDM_21_species, pars = "Lambda")
#L.sam <- as.matrix(fit.beta.JSDM.n500, pars = "Lambda")

L.post_means <- colMeans(L.sam)
# reorder as matrix (n_latent_factors x n_species)
L.post_means <- matrix(L.post_means, ncol = ncol(Y), byrow = FALSE)

plot(L.post_means[1,],L.post_means[2,], pch = 19, xlab = "loading 1", ylab="loading 2")
text(L.post_means[1,],L.post_means[2,], colnames(Y), cex = 0.8, pos = 3)
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# one can also calculate the residual covariance matrix
cor.mat <- cov2cor(t(L.post_means) %*% L.post_means)
colnames(cor.mat) <- rownames(cor.mat) <- colnames(Y)
cor.mat


### examine the convergence

fit.beta.JSDM.21_species <- readRDS("models/multivariate/n_500/M1/JSDM_21_species.RDS")

Rhats <- summary(fit.beta.JSDM.21_species)$summary[,"Rhat"]
sort(Rhats,decreasing = TRUE) # there are some non-converging chains

# examine the trace plots
stan_trace(fit.beta.JSDM.21_species, pars = paste0("Lambda[2,",6:21,"]"))

### check the convergence in terms of covariance matrix S = L^T %*% L
post.samples <- extract(fit.beta.JSDM.21_species)
post.samples.array <- as.array(fit.beta.JSDM.21_species)

n_species <- ncol(Y)
sp_names <- colnames(Y)
n_factors <- 2
n_samples <- dim(post.samples.array)[1]
n_chains <- dim(post.samples.array)[2]
R_hats <- c()
res_list <- list("mean(R)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)),
                 "Pr(R>0)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)),
                 "Pr(R<0)" = matrix(0,nrow=n_species,ncol=n_species,dimnames = list(sp_names,sp_names)))

names_of_rows <- c()
par(mfrow = c(5,5),
    mar = c(2,2,2,0))
for(i in 1:(n_species-1)) {
  for (j in (i+1):n_species) {
    Sigma_ij <- matrix(0,nrow=n_samples,ncol=n_chains)
    Sigma_ii <- matrix(0,nrow=n_samples,ncol=n_chains)
    Sigma_jj <- matrix(0,nrow=n_samples,ncol=n_chains)
    for (k in 1:n_factors) {
      L_ki <- post.samples.array[,,paste0("Lambda[",k,",",i,"]")]
      L_kj <- post.samples.array[,,paste0("Lambda[",k,",",j,"]")]
      Sigma_ij <- Sigma_ij + L_ki*L_kj
      Sigma_ii <- Sigma_ii + L_ki^2
      Sigma_jj <- Sigma_jj + L_kj^2
    }
    # from covariance to correlation
    R_ij <- Sigma_ij / (sqrt(Sigma_ii)*sqrt(Sigma_jj))
    
    # calculate R-hat
    Rh <- Rhat(R_ij)
    R_hats[paste0("S_",i,j)] <- Rh
    
    # save posterior mean, Pr(R_ij > 0) and Pr(R_ij < 0)
    res_list[[1]][i,j] <- mean(R_ij)
    res_list[[2]][i,j] <- mean(R_ij > 0)
    res_list[[3]][i,j] <- mean(R_ij < 0)

    # plot chains
    matplot(Sigma_ij, type = "l", lty = 1, lwd = 2, col = 1:n_chains)
    legend("bottomleft", legend = paste0("Rhat: ", round(Rh,2)), bty = "n", col = "red")
  }
}

### try to visualize correlation matrix

# 1) just posterior means
par(mfrow = c(1,1))
corrplot(res_list$`mean(R`, type = "upper", order = "original")

# 2) different colors for 1) positive correlation 2) negative correlation 3) no correlation (overlaps with 0)

R_mat_categorical <- matrix(0,nrow=n_species,ncol=n_species,dimnames=list(sp_names,sp_names))
R_mat_categorical[which(res_list$`Pr(R>0` > 0.95)] <- 1
R_mat_categorical[which(res_list$`Pr(R<0` > 0.95)] <- -1
corrplot(R_mat_categorical, type = "upper", order = "original")

