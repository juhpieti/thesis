#### spatial versions of model 3 and 4
#### fit for two species
#### compare the result with non-spatial

# load in pacakges
library(terra)
library(loo)
library(ggplot2)

# load in stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the training data
#load("data/estonia_sub/estonia_sub_df.RData")
load("data/estonia_new/train_2020_2021.Rdata")

# load in the helper functions
source("codes/helpers.R")

estonia_sub <- df_sub

# in the data there are three species with less than 5 observations, I will delete them
estonia_sub <- estonia_sub[,!(colnames(estonia_sub) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]


### visualize the data spatially
europe.vect <- vect("europe_coastline_shapefile/Europe_coastline_poly.shp")
europe.vect <- project(europe.vect, "EPSG:3067") #reproject to TM35FIN

# crop for some region of interest
ext_vec <- c(0,566613.8,6000000,6650000)
europe.vect.cut <- crop(europe.vect,ext_vec)

estonia_sub.vect <- vect(estonia_sub, geom = c("x","y"), crs = "EPSG:3067")

par(mfrow = c(1,1))
plot(estonia_sub.vect, cex = 1, col = "red", main = "observations")
plot(europe.vect.cut, add = TRUE, col = "lightgrey")
plot(estonia_sub.vect, cex = 1, col = "red", main = "observations", add = TRUE)

# ### for spatial random effect we want an even grid to the area except the coast
# # parameters of the grid
# #ext_vec_grid <- ext(c(120000,560000,6300000,6640000))
# ext_vec_grid <- ext(c(140000,560000,6320000,6640000))
# length_grid_cell <- 20000 #20km x 20km grid
# ncols <- ( ext_vec_grid[2] - ext_vec_grid[1] ) / length_grid_cell
# nrows <- ( ext_vec_grid[4] - ext_vec_grid[3] ) / length_grid_cell
# 
# grid_raster <- rast(ext_vec_grid,ncols=ncols,nrows=nrows, crs = "EPSG:3067")
# #values(grid_raster) <- sample(1:ncell(grid_raster))
# 
# # erase the land areas
# grid_vector <- as.polygons(grid_raster)
# spatial_grid <- erase(grid_vector, europe.vect.cut)
# 
# # plot the grid
# par(mfrow = c(1,1))
# plot(spatial_grid, main = "prediction grid")
# plot(estonia_sub.vect, add = TRUE, col = "red")
# plot(europe.vect.cut, add = TRUE, col = "lightgrey")
# 
# # save the spatial random effect grid
# writeVector(spatial_grid, filename = "data/estonia_new/spatial_random_effect_grid.shp", overwrite = TRUE)

### take centre points and create space matrix (s1,s2) with coordinates
### put it in the model such that they have a gaussian process
### in likelihood part, take spatial errors using PZ type of thing that picks the correct random effect

### load in spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]

dim(grid_centers.df) # there are m = 191 grid cells

### create a P matrix (nxm) matrix that tells to which grid cell each observation belongs to
nearest_grid_center <- nearest(estonia_sub.vect, grid_centers)
nearest_grid_center.df <- as.data.frame(nearest_grid_center)

# n-length vector indicating the ID of the nearest grid center
nearest_grid_center.vec <- nearest_grid_center.df$to_id

# take the indexes of grid cells that has observations in them
observed_grid_cells <- unique(nearest_grid_center.vec)
observed_grid_cells.df <- grid_centers.df[observed_grid_cells,c("x","y")]

# create the P matrix
P <- matrix(0,ncol=length(observed_grid_cells),nrow=nrow(estonia_sub))
colnames(P) <- rownames(observed_grid_cells.df)
for (i in 1:nrow(estonia_sub)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1 
}

### MODELING ###

X <- estonia_sub[,11:19]
### scale the covariates
X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

### Amphibalanus & Chara Aspera

amphi <- estonia_sub[,"Amphibalanus improvisus"]
amphi.01 <- amphi/100 # beta regression takes values from [0,1]

chara <- estonia_sub[,"Chara aspera"]
chara.01 <- chara/100


loo_amphi.spat <- c()

### put the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

### run models
### spatial (without zero-inflation)
amphi_dat.beta_spat <- list(N = nrow(X.sec_ord),
                       n_var = ncol(X.sec_ord),
                       n_obs_grid = ncol(P),
                       y = amphi.01,
                       X = X.sec_ord,
                       s = observed_grid_cells.df,
                       a = 1,
                       P = P)



mod_amphi.beta_spat <- stan("stan_files/left_censored_beta_regression_spatial.stan", data = amphi_dat.beta_spat, chains = 4, iter = 200, seed = 42,
                            pars = c("mu","logneg_beta_2"), include = FALSE)


saveRDS(mod_amphi.beta_spat, "models/demo/M3/amphi.rds")
saveRDS(mod_chara.beta_spat, "models/demo/M3/chara.rds")



#save loo
loo_amphi.spat <- c(loo_amphi.spat,loo(mod_amphi.beta_spat)$elpd_loo)

### posterior predictive checks
# histogram of coverages in original data
par(mfrow = c(4,4),
    mar = c(2,4,2,0))
hist(amphi.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_amphi.beta_spat, pars = c("alpha"))
beta.sam <- as.matrix(mod_amphi.beta_spat, pars = c("beta"))
rho.sam <- as.matrix(mod_amphi.beta_spat, pars = c("rho"))
phi.sam <- as.matrix(mod_amphi.beta_spat, pars = c("phi"))

a <- 1

# make n_rep replications of new data, compare to observed coverages
n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)
for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,] + P[i,] %*% phi.sam[idx,])
    rho_i <- rho.sam[idx,]
    V_i <- rbeta(1,mu_i*rho_i,(1-mu_i)*rho_i)
    y_rep_i <- max(0,(a+1)*V_i - a)
    y_rep <- c(y_rep, y_rep_i)
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_beta_regression(mod_amphi.beta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"expectation")
plot_responses_beta_regression(mod_chara.beta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"expectation")

### plot prob of zero
plot_responses_beta_regression(mod_amphi.beta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"probzero")

### plot distribution of y in average conditions
average_distribution_beta_regression(mod_amphi.beta_spat,1,"amphi")

### plot the spatial random effects!
s_obs <- observed_grid_cells.df
s_pred <- grid_centers.df[,c("x","y")]/1000 #as kilometers
grid.vect <- spatial_grid
shoreline.vect <- europe.vect.cut
obs.vect <- estonia_sub.vect

# plot the spatial effects
plot_spatial_effects_beta(mod_amphi.beta_spat,s_obs,s_pred,grid.vect,obs.vect,shoreline.vect)
  
#save the loo values
load("results/demo/loo_table.Rdata")
loo_table <- cbind(loo_table, "Beta spat." = c(loo_amphi.spat,rep(0,3)))
loo_table


### Spatial ZI-Beta Regression

mod_amphi.ZIBeta_spat <- stan("stan_files/zero_inflated_left_censored_beta_regression_spatial.stan", data = amphi_dat.beta_spat, chains = 4, iter = 200, seed = 42,
                              pars = c("mu","prob_suit","logneg_beta_2","logneg_beta_pi_2"), include = FALSE)
saveRDS(mod_amphi.ZIBeta_spat, "models/demo/M4/amphi.rds")


mod_chara.ZIBeta_spat <- stan("stan_files/zero_inflated_left_censored_beta_regression_spatial.stan", data = chara_dat.beta_spat, chains = 4, iter = 200, seed = 42,
                              pars = c("mu","prob_suit","logneg_beta_2","logneg_beta_pi_2"), include = FALSE)


loo_amphi.spat <- c(loo_amphi.spat,loo(mod_amphi.ZIBeta_spat)$elpd_loo)

#posterior predictive
par(mfrow = c(4,4),
    mar = c(2,4,2,0))

hist(amphi.01, breaks = 10, xlim = c(0,1), main = "obs", ylim = c(0,500))

alpha.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("alpha"))
alpha_pi.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("alpha_pi"))
beta.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("beta"))
beta_pi.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("beta_pi"))
rho.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("rho"))
phi_mu.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("phi_mu"))
phi_pi.sam <- as.matrix(mod_amphi.ZIBeta_spat, pars = c("phi_pi"))

a <- 1

n_rep <- 15
rep_idx <- sample(1:nrow(beta.sam),n_rep,replace = FALSE)

for (i in 1:n_rep) {
  idx <- rep_idx[i]
  y_rep <- c()
  for (j in 1:nrow(X.sec_ord)) {
    pi_i <- inv_logit(alpha_pi.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta_pi.sam[idx,] + P[i,] %*% phi_pi.sam[idx,])
    Z_i <- rbinom(1,1,pi_i)
    if (Z_i == 0) {
      y_rep <- c(y_rep,0)
    } else {
      mu_i <- inv_logit(alpha.sam[idx,] + as.numeric(X.sec_ord[j, ]) %*% beta.sam[idx,] + P[i,] %*% phi_mu.sam[idx,])
      rho_i <- rho.sam[idx,]
      V_i <- rbeta(1,mu_i*rho_i,(1-mu_i)*rho_i)
      y_rep_i <- max(0,(a+1)*V_i - a)
      y_rep <- c(y_rep, y_rep_i)
    }
  }
  hist(y_rep, breaks = 10, xlim = c(0,1), main = paste0("rep",i), ylim = c(0,500))
}

### plot response curves
plot_responses_ZI_beta_regression(mod_amphi.ZIBeta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"expectation")

### plot prob of zero
plot_responses_ZI_beta_regression(mod_amphi.ZIBeta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,1,3,4,1,"probzero")

### plot distribution of y in average conditions
average_distribution_ZI_beta_regression(mod_amphi.ZIBeta_spat,1,"amphi")

### plot separate curves
plot_separate_responses_ZI_beta_regression(mod_amphi.ZIBeta_spat,X.sec_ord,X,TRUE,TRUE,0,60,0,100,3,4,1)


### plot the spatial random effects!!! (two of them)

s_obs <- observed_grid_cells.df
s_pred <- grid_centers.df[,c("x","y")]/1000 #as kilometers
grid.vect <- spatial_grid
shoreline.vect <- europe.vect.cut
obs.vect <- estonia_sub.vect

plot_spatial_effects_ZIBeta(mod_amphi.ZIBeta_spat,s_obs,s_pred,grid.vect,obs.vect,shoreline.vect)

#save the loo values
loo_table <- cbind(loo_table, "ZI-Beta spat." = rep(0,4))

save(loo_table, file = "results/demo/loo_table.Rdata")
