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
train <- train_n100

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
Y <- Y[,colSums(Y>0)>10] #take only species with > 10 appearances

cor.mat <- cor(Y)
corrplot(cor.mat, type = "upper", order = "hclust", method = "number")

# take interesting, somewhat correlated species
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y <- Y[ ,sp_names]

### SPATIAL 

### load in spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]

dim(grid_centers.df) # there are m = 182 grid cells

### create a P matrix (nxm) matrix that tells to which grid cell each observation belongs to
estonia_sub.vect <- vect(train, geom = c("x","y"), crs = "EPSG:3067")

nearest_grid_center <- nearest(estonia_sub.vect, grid_centers)
nearest_grid_center.df <- as.data.frame(nearest_grid_center)

# n-length vector indicating the ID of the nearest grid center
nearest_grid_center.vec <- nearest_grid_center.df$to_id

# take the indexes of grid cells that has observations in them
observed_grid_cells <- unique(nearest_grid_center.vec)
observed_grid_cells.df <- grid_centers.df[observed_grid_cells,c("x","y")]

# create the P matrix
P <- matrix(0,ncol=length(observed_grid_cells),nrow=nrow(train))
colnames(P) <- rownames(observed_grid_cells.df)
for (i in 1:nrow(X)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1 
}

### put the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000


### MODELING PART
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms

Y.scaled <- Y/100 #scale Y between 0 and 1

n_chains <- 4
n_iter <- 2000

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

fit.beta.JSDM.spatial <- stan("stan_files/multivariate/left_censored_beta_regression_spatial_multivariate.stan",
                              data = data.spat.list, chains = n_chains, iter = n_iter, seed = 42,
                              pars = c("Mu"), include = FALSE)

### SAVE THE MODEL FIT
subfolder <- paste0("n_",nrow(Y.scaled),"/",ncol(X),"_covariates/")
saveRDS(fit.beta.JSDM.spatial, file = paste0("models/multivariate/",subfolder,"JSDM_spatial_test.RDS"))

### EXAMINE THE MODEL FIT
fit.beta.JSDM.spatial <- readRDS("models/multivariate/JSDM_spatial_test.RDS")
loo(fit.beta.JSDM.spatial)

### SITE LOADINGS (color with some missing variables?)
Z.sam <- as.matrix(fit.beta.JSDM.spatial, pars = "phi")
Z.post_means <- colMeans(Z.sam)
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


### SPECIES LOADINGS
L.sam <- as.matrix(fit.beta.JSDM.spatial, pars = "Lambda")
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
