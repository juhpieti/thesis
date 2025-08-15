#######################################################################################
### SCRIPT TO RUN MULTIVARIATE SPATIAL ZERO-INFLATED LEFT-CENSORED BETA REGRESSION ####
#######################################################################################

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

### PREPARE THE DATA FOR SPATIAL MODELING

### load in the coarse spatial random effect grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

# prepare a matrix of grid centers and their locations (TM35FIN coordinates)
grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]
dim(grid_centers.df) # there are m = 191 grid cells

# turn the training data into terra-vector
train.vect <- vect(train, geom = c("x","y"), crs = "EPSG:3067")

# find the spatial grid cell for each sampling location
nearest_grid_center <- nearest(train.vect, grid_centers)
nearest_grid_center.df <- as.data.frame(nearest_grid_center)

# n-length vector indicating the ID of the nearest grid center
nearest_grid_center.vec <- nearest_grid_center.df$to_id

# take the indexes of grid cells that have observations in them
observed_grid_cells <- unique(nearest_grid_center.vec)
observed_grid_cells.df <- grid_centers.df[observed_grid_cells,c("x","y")]

# create matrix P (N x n_grid_cells), where i:th row indicates the spatial grid cell that i:th sampling point is located in
# so each row of matrix P sums up to 1
P <- matrix(0,ncol=length(observed_grid_cells),nrow=nrow(train))
colnames(P) <- rownames(observed_grid_cells.df)
for (i in 1:nrow(estonia_sub)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

### turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

### MODELING PART ###
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms
Y.scaled <- Y/100 #scale y to (0,1) for beta regression

# stan input parameters
n_chains <- 4
n_iter <- 2000

# prepare the data for stan
data.spat.list <- list(N = nrow(Y.scaled),
                       n_var = ncol(X.sec_ord),
                       n_obs_grid = nrow(observed_grid_cells.df),
                       J = ncol(Y.scaled),
                       n_f = 2,
                       Y = Y.scaled,
                       X = X.sec_ord,
                       s = observed_grid_cells.df,
                       a = 1,
                       P = P)

fit.ZIbeta.JSDM.spatial <- stan("stan_files/multivariate/zero_inflated_left_censored_beta_regression_spatial_multivariate.stan",
                              data = data.spat.list, chains = n_chains, iter = n_iter, seed = 42,
                              pars = c("Mu","Pi"), include = FALSE)

### SAVE THE MODEL FIT
subfolder <- paste0("n_",nrow(Y.scaled),"/",ncol(X),"_covariates/")
saveRDS(fit.ZIbeta.JSDM.spatial, file = paste0("models/multivariate/",subfolder,"JSDM_spatial_ZI_test.RDS"))

### EXAMINE THE MODEL FIT
fit.ZIbeta.JSDM.spatial <- readRDS("models/multivariate/JSDM_spatial_ZI_test.RDS")
loo(fit.beta.JSDM.spatial)

### VISUALIZE SITE LOADINGS (for mean of beta distribution)
Z.sam <- as.matrix(fit.ZIbeta.JSDM.spatial, pars = "phi")
Z.post_means <- colMeans(Z.sam)
# reorder as matrix (N x n_latent_factors)
Z.post_means <- matrix(Z.post_means, nrow = nrow(observed_grid_cells.df), byrow = FALSE)

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

### VISUALIZE SITE LOADINGS (for prob. of suitability)
Z_pi.sam <- as.matrix(fit.ZIbeta.JSDM.spatial, pars = "phi_pi")
Z_pi.post_means <- colMeans(Z_pi.sam)
Z_pi.post_means <- matrix(Z_pi.post_means, nrow = nrow(observed_grid_cells.df), byrow = FALSE)

plot(Z_pi.post_means[,1],Z_pi.post_means[,2], pch = 19, xlab = "latent factor 1", ylab = "latent factor 2")
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# Color points wrt. covariates not used in modeling
df_pi <- as.data.frame(cbind(Z_pi.post_means,train_n100[,c("no3_bottom","po4_bottom","current_bottom","chl_bottom")]))
colnames(df_pi)[1:2] <- c("lf1","lf2")

plots <- list()
for (var in colnames(df_pi)[3:6]) {
  plots[[var]] <- ggplot(data=df_pi,aes(x = lf1, y = lf2)) + 
    geom_point(aes(colour=.data[[var]]))
}
gridExtra::grid.arrange(grobs = plots, ncol = 2)


### SPECIES LOADINGS (mean of beta distribution)
L.sam <- as.matrix(fit.ZIbeta.JSDM.spatial, pars = "Lambda")
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

### SPECIES LOADINGS (prob. of suitability)
L_pi.sam <- as.matrix(fit.ZIbeta.JSDM.spatial, pars = "Lambda_pi")
L_pi.post_means <- colMeans(L_pi.sam)
# reorder as matrix (n_latent_factors x n_species)
L_pi.post_means <- matrix(L_pi.post_means, ncol = ncol(Y), byrow = FALSE)

plot(L_pi.post_means[1,],L_pi.post_means[2,], pch = 19, xlab = "loading 1", ylab="loading 2")
text(L_pi.post_means[1,],L_pi.post_means[2,], colnames(Y), cex = 0.8, pos = 3)
abline(v = 0, col = "red", lty = 2)
abline(h = 0, col = "red", lty = 2)

# one can also calculate the residual covariance matrix
cor.mat <- cov2cor(t(L_pi.post_means) %*% L_pi.post_means)
colnames(cor.mat) <- rownames(cor.mat) <- colnames(Y)
cor.mat
