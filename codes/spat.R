### MODEL 1: left-censored beta-regression with spatial random effects
### this script
# 1) fits
# 2) plots responses
# 3) plots predicted maps
# 4) ?

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
load("data/estonia_new/train_2020_2021.Rdata")
colnames(df_sub)
colSums(df_sub[,20:38] > 0)

# remove species that are too rare (under 5 observations)
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

# PREPARE THE COORDINATES

# load in the predictive grid
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")

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
for (i in 1:nrow(estonia_sub)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1 
}

### put the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

# loop over the species, save the models
sp_names <- colnames(train)[20:35]
n_chains <- 4
n_iter <- 50
for (sp_name in sp_names[1:1]) {
  y <- train[,sp_name]
  y.01 <- y/100
  
  dat.beta_spat <- list(N = nrow(X.sec_ord),
                              n_var = ncol(X.sec_ord),
                              n_obs_grid = ncol(P),
                              y = y.01,
                              X = X.sec_ord,
                              s = observed_grid_cells.df,
                              a = 1,
                              P = P)
  
  
  
  mod.beta_spat <- stan("stan_files/left_censored_beta_regression_spatial.stan", data = dat.beta_spat, chains = n_chains, iter = n_iter, seed = 42,
                              pars = c("mu","logneg_beta_2"), include = FALSE)
  
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  f_name <- paste0("models/M3/",sp_name_modified,".rds")
  
  saveRDS(mod.beta_spat, f_name)
}



#load in the models
mod_amphi.spat <- readRDS("models/M3/Amphibalanus_improvisus.rds")
mod_amphi.spat <- readRDS("models/demo/M3/amphi.rds")
mod_chara.spat <- readRDS("models/demo/M3/chara.rds")

pred_list_m3 <- predict_spatial_beta_regression(mod_amphi.spat,pred_grid_1km_2021_july_df[,2:10],X,
                                                     pred_grid_1km_2021_july_df[,c("x","y")],
                                                     grid_centers.df/1000, observed_grid_cells.df,10,1)

pred_list_m3 <- predict_spatial_beta_regression(mod_amphi.spat,
                                                cbind(pred_grid_1km_2021_july_df[,2:10],pred_grid_1km_2021_july_df$depth / pred_grid_1km_2021_july_df$zsd),
                                                X,
                                                pred_grid_1km_2021_july_df[,c("x","y")],
                                                grid_centers.df/1000, observed_grid_cells.df,10,1)



locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m3$EY_sam)
vars <- apply(pred_list_m3$EY_sam,2,var)
#vals <- colMeans(pred_list_spatial$y_sam)
vals.scaled <- vals/max(vals)
test_df <- as.data.frame(cbind(locs, vals, vars, vals.scaled))

### try the acumulative thing
head(test_df)
test_df_ordered <- test_df[order(-test_df$vals), ]
test_df_ordered$hotspot <- 1

# find when 90% of the expected value is reached
tot_sum <- sum(test_df_ordered$vals)
idx <- which(cumsum(test_df_ordered$vals) > 0.7*tot_sum)[1]

# set all outside to be non-hotspot
test_df_ordered[idx:nrow(test_df_ordered),"hotspot"] <- 0


test_vect <- vect(test_df_ordered, geom=c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r <- rasterize(test_vect, test_rast, field = "vals")
r.vars <- rasterize(test_vect, test_rast, field = "vars")
r.hotspot <- rasterize(test_vect, test_rast, field = "hotspot")
r.scaled <- rasterize(test_vect, test_rast, field = "vals.scaled") 

plot(r, colNA = "lightgrey", main = "M3 (w/o ZI,  w/ RE)",plg = list(title = "E[coverage]"),
     xlab = "Easting (m)", ylab = "Northing (m)")
#plot(r.scaled, colNA = "lightgrey", main = "Amphi rel. exp. coverage: M3 (w/o ZI,  w/ RE)")
plot(r.hotspot, colNA = "lightgrey", main = "M3 (w/o ZI, w/ RE)", col = c("red","blue"),
     plg = list(title = "hotspot", legend = c("no","yes")),
     xlab = "Easting (m)", ylab = "Northing (m)")

# plot variance
plot(r.vars,colNA="lightgrey", main = "Amphi variance: M3 (w/o ZI, w/ RE)")


### plot the random effects
locs <- pred_grid_1km_2021_july_df[,c("x","y")]
vals <- colMeans(pred_list_m3$phi_pred_sam)
test_df <- as.data.frame(cbind(locs, vals))

test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.phi <- rasterize(test_vect, test_rast, field = "vals")

par(mfrow = c(1,1))
plot(r.phi, colNA = "lightgrey", main = "Amhpi spatial REs: M3 (w/o ZI, w/ RE)")
plot(r.phi, colNA = "lightgrey", main = "Chara spatial REs: M3 (w/o ZI, w/ RE)")


### difference with M1
vals <- colMeans(pred_list_spatial$EY_sam) - colMeans(pred_list$EY_sam)
#vals <- colMeans(pred_list_spatial$y_sam) - colMeans(pred_list$y_sam)
test_df <- as.data.frame(cbind(locs, vals))
test_vect <- vect(test_df, geom = c("x","y"), crs = "EPSG:3067")
test_rast <- rast(ext = ext(predictive_grid), res = 1000, crs = "EPSG:3067")
r.diff <- rasterize(test_vect, test_rast, field = "vals")

par(mfrow = c(1,1))
plot(r.diff, colNA = "lightgrey", main = "Amphi M3-M1 expected coverage")
plot(r.diff, colNA = "lightgrey", main = "Chara M3-M1 expected coverage")


