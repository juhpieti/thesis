#########################################################
### SCRIPT TO RUN 2-STAGE JSDM FOR PERCENT COVER DATA ###
#########################################################

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

### load in the dataset(s)

load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
dim(train_n100_all_species)
train <- train_n100_all_species

load("data/estonia_new/train/train_2020_2021_all_species_n250.Rdata")
dim(train_n250_all_species)
train <- train_n250_all_species

load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
dim(train_n500_all_species)
train <- train_n500_all_species

### FULL DATA (n = 3221)
load("data/estonia_new/train/train_2020_2021_full_all_species.Rdata")
dim(df_all_species)
colnames(df_all_species)
sort(colSums(df_all_species[,20:71] > 0), decreasing = TRUE) # still quite many species are not observed at all
sum(colSums(df_all_species[,20:71] > 0) > 0) # so 27/51 species are observed in these 3221 sites
train <- df_all_species

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
for (i in 1:nrow(train)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

### turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

### VISUALIZE THE DATA
## compare how total mass compares with 4 species to 21 species
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
Y_4 <- train[,sp_names]

colnames(train)
Y <- train[,20:71]
idx_positive <- colSums(Y > 0) > 2 # take only species that have some observations
sum(idx_positive) # it's 21 of those
Y_21 <- Y[,idx_positive]
Y <- Y_21[,!colnames(Y_21) == "Ranunculus peltatus subsp_ Baudotii"] # drop 1 species for convenient J = 20

### see the distributions both on normal scale and log-scale
total_mass_4sp <- rowSums(Y_4)
total_mass_21sp <- rowSums(Y_21)

par(mfrow = c(2,2))
hist(total_mass_4sp, breaks=25)
hist(total_mass_21sp, breaks = 25)
hist(log(total_mass_4sp), breaks = 25)
hist(log(total_mass_21sp), breaks = 25)

### examine how the mass behaves wrt. covariates
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) #approximate the light level at the bottom
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

par(mfrow = c(3,3),
    mar = c(4,4,2,0))
for (i in 1:ncol(X)) {
  plot(X[,i], total_mass_21sp, pch = 18, col = "red", xlab = colnames(X)[i], ylab = "total mass")
  points(X[,i], total_mass_4sp, pch = 18, col = "blue")
  if (i == 1) (legend("topright", legend = c("J=4","J=21"), col = c("blue","red"), pch = c(18,18)))
}

par(mfrow = c(3,3),
    mar = c(4,4,2,0))
for (i in 1:ncol(X)) {
  plot(X[,i], log(total_mass_21sp), pch = 18, col = "red", xlab = colnames(X)[i], ylab = "log(total mass)")
  points(X[,i], log(total_mass_4sp), pch = 18, col = "blue")
  if (i == 1) (legend("topright", legend = c("J=4","J=21"), col = c("blue","red"), pch = c(18,18)))
}

### fit the model
X.scaled <- scale_covariates(X) #scale the covariates
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled)) #add second order terms

# prepare data for stan
data.spat.list <- list(N = nrow(Y),
                       n_var = ncol(X.sec_ord),
                       n_obs_grid = nrow(observed_grid_cells.df),
                       J = ncol(Y),
                       n_f = 2,
                       Y = Y,
                       y_sum = rowSums(Y),
                       X = X.sec_ord,
                       s = observed_grid_cells.df,
                       P = P)

data.spat.list.nf3 <- list(N = nrow(Y),
                       n_var = ncol(X.sec_ord),
                       n_obs_grid = nrow(observed_grid_cells.df),
                       J = ncol(Y),
                       n_f = 3,
                       Y = Y,
                       y_sum = rowSums(Y),
                       X = X.sec_ord,
                       s = observed_grid_cells.df,
                       P = P)

# stan input parameters
n_chains <- 4
n_iter <- 2000

# fit the model
# fit.2stage.poisson.DirMult.spat <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_poisson_spatial.stan",
#                                    data = data.spat.list, chains = n_chains, iter = n_iter, seed = 42,
#                                    pars = c("Mu_Dir", "mu_M"), include = FALSE)

fit.2stage.NegBin.DirMult.spat <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_NegBin_spatial.stan",
                                        data = data.spat.list, chains = n_chains, iter = n_iter, seed = 42,
                                        pars = c("Mu_Dir", "mu_M", "z_M","Z_spat"), include = FALSE)

# fit.2stage.NegBin.DirMult.spat.hier_priors <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_NegBin_spatial_hier_priors.stan",
#                                        data = data.spat.list, chains = n_chains, iter = n_iter, seed = 42,
#                                        pars = c("Mu_Dir", "mu_M"), include = FALSE)


fit.2stage.NegBin.DirMult.spat.nf3 <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_NegBin_spatial.stan",
                                       data = data.spat.list.nf3, chains = n_chains, iter = n_iter, seed = 42,
                                       pars = c("Mu_Dir", "mu_M", "z_M", "Z_spat"), include = FALSE)

# fit.2stage.poisson.DirMult.hier_priors.spat <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_poisson_hier_priors.stan",
#                                                data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                                                pars = c("Mu_Dir", "mu_M"), include = FALSE)
# 
# 
# fit.2stage.NegBin.DirMult.spat <- stan("stan_files/multivariate/2stage_JSDM_DirMultinomial_negbin.stan",
#                                   data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                                   pars = c("Mu_Dir", "mu_M"), include = FALSE)



# fit.2stage <- stan("stan_files/multivariate/2stage_JSDM.stan",
#                       data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                       pars = c("mu"), include = FALSE)
# 
# fit.2stage.gamma <- stan("stan_files/multivariate/2stage_JSDM_gamma.stan",
#                    data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                    pars = c("mu"), include = FALSE)
# 
# fit.2stage.beta <- stan("stan_files/multivariate/2stage_JSDM_scaled_beta.stan",
#                          data = data.list, chains = n_chains, iter = n_iter, seed = 42,
#                          pars = c("mu"), include = FALSE)

# save the models
subfolder <- paste0("n_",nrow(Y),"/")
saveRDS(fit.2stage.NegBin.DirMult.spat, paste0("models/two_stage/",subfolder,"M3/JSDM_NegBin_DirMult_spatial.RDS"))
#saveRDS(fit.2stage.NegBin.DirMult.spat.hier_priors, paste0("models/two_stage/",subfolder,"M3/JSDM_NegBin_DirMult_spatial_hier_priors.RDS"))