############################################################
### SCRIPT TO ANALYSE RESULTS FROM FITTED 2-STAGE JSDM   ###
### INCLUDES PLOTTING MAPS, DRAWING RESPONSE CURVES      ###
############################################################

### models fitted with separate R-scripts (joint_beta_regression.R, joint_zi_beta_regression.R, joint_beta_regression_spatial.R, joint_zi_beta_regression_spatial.R)
### this script loads fitted models and draws different maps and curves, as well as produces tables
library(terra)
library(loo)

# load in helper functions
source("codes/helpers.R")

# load in the training data
# load("data/estonia_new/train/train_2020_2021_n100.Rdata")
# train <- train_n100

load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
train <- train_n100_all_species

# load("data/estonia_new/train/train_2020_2021_n500.Rdata")
# train <- train_n500

colnames(train)
colSums(train[,20:38] > 0)

# remove species that are too rare (under 5 observations)
#train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
#X$depth_to_secchi <- X$depth / X$zsd # add secchi/depth for a variable representing seafloor light level
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

Y <- train[,20:71]
idx_positive <- colSums(Y > 0) > 0 # take only species that have some observations
sum(idx_positive) # it's 21 of those
Y_21 <- Y[,idx_positive]
Y <- Y_21[,!colnames(Y_21) == "Ranunculus peltatus subsp_ Baudotii"] # drop 1 species for convenient J = 20

# load in the predictive grid
predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

# add relative light level
pred_grid_1km_2021_july_df$light_bottom <- exp(-1.7*pred_grid_1km_2021_july_df$depth/pred_grid_1km_2021_july_df$zsd)

### load in spatial grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]

dim(grid_centers.df) # there are m = 191 grid cells

### find the observed grid cells
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
for (i in 1:nrow(train)) {
  P[i,as.character(nearest_grid_center.vec[i])] <- 1
}

### put the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

### prepare functions

plot_random_effects_2stage <- function(pred_list,locs,pred.grid.vect,title_chr="") {
  ### type: "mu" or "pi"
  
  # how many latent factors are there in the model?
  n_factors <- dim(pred_list$Z_mu_sam)[3]
  
  for (k in 1:n_factors) {
    vals <- colMeans(pred_list$Z_mu_sam[,,k])
    pred_df <- as.data.frame(cbind(locs,vals))
    pred_vect <- vect(pred_df, geom = c("x","y"), crs = "EPSG:3067")
    pred_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
    r <- rasterize(pred_vect, pred_rast, field = "vals")
    
    plot(r, colNA = "lightgrey", main = title_chr, plg = list(title = paste0("Z",k)), xlab = "Easting (m)", ylab = "Northing (m)")
  }
  
  # plot also the random effects for total cover
  vals <- colMeans(pred_list$phi_M_sam)
  pred_df <- as.data.frame(cbind(locs,vals))
  pred_vect <- vect(pred_df, geom = c("x","y"), crs = "EPSG:3067")
  pred_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(pred_vect, pred_rast, field = "vals")
  
  plot(r, colNA = "lightgrey", main = title_chr, plg = list(title = "phi_Ytot"), xlab = "Easting (m)", ylab = "Northing (m)")
}



### DRAW MAPS ###
sp_names <- names(JSDM_predictions)




### 2) JSDM

# load model (20 species)
fit.JSDM.20species <- readRDS("models/multivariate/n_100/M1/JSDM_20species_hier_priors.RDS")

# predict (J=4)
sp_names <- colnames(Y)

set.seed(42)
JSDM_predictions_20species <- predict_beta_regression_JSDM(fit.JSDM.20species,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,100,1,FALSE,1000,0,FALSE,0)
plot_map_JSDM(JSDM_predictions_20species,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = TRUE)



# predict (J=21)
sp_names <- colnames(train_n100_all_species[,20:71])
sp_names <- sp_names[which(colSums(train_n100_all_species[,20:71] > 0) > 0)]

set.seed(42)
JSDM_predictions <- predict_beta_regression_JSDM(fit.JSDM.21species,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,100,1,FALSE,1000,0,FALSE,0)
set.seed(42)
JSDM_predictions.hier.prior <- predict_beta_regression_JSDM(fit.JSDM.21species.hier.prior,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,100,1,FALSE,1000,0,FALSE,0)

# plot
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = FALSE)
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = TRUE)
plot_map_JSDM(JSDM_predictions.hier.prior,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = TRUE)


# plot also hotspots
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.8,sp_names,summary_maps = FALSE)
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.5,sp_names,summary_maps = TRUE)


##### 2) ZERO-INFLATED LEFT-CENSORED BETA REGRESSION
plot_map_JSDM(JSDM.ZI_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)

##### 3) SPATIAL LEFT-CENSORED BETA REGRESSION
# gather the predictions as lists
stacked_SDM_spat_predictions <- list()
for(sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  pred <- predict_spatial_beta_regression(mod,pred_grid_1km_2021_july_df[,colnames(X)],X,pred_grid_1km_2021_july_df[,c("x","y")],grid_centers.df/1000, observed_grid_cells.df,400,1,FALSE,1000)
  stacked_SDM_spat_predictions[[sp_name]] <- pred
}


# plot with JSDM function
plot_map_JSDM(stacked_SDM_spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = FALSE)
plot_map_JSDM(stacked_SDM_spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = TRUE)

# plot also hotspots
plot_map_JSDM(stacked_SDM_spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.8,sp_names,summary_maps = FALSE)
plot_map_JSDM(stacked_SDM_spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.5,sp_names,summary_maps = TRUE)

# plot JSDMs
plot_map_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = FALSE)
plot_map_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,summary_maps = TRUE)

# plot also hotspots
plot_map_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.8,sp_names,summary_maps = FALSE)
plot_map_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.5,sp_names,summary_maps = TRUE)


# plot spatially correlated latent factors
plot_random_effects_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","latent factors (mu)")

### loop over 4 species, predict and see the random effects
plot_random_effects(SDM_spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","RE (mu)")


##### 4) ZERO-INFLATED LEFT-CENSORED BETA REGRESSION
plot_map_JSDM(JSDM.spat_ZI_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)



### Load in models

### Plot coverages
png(paste0("plots/final_results/",subfolder,"/coverage_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
par(mfrow = c(2,2))
plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE")
plot_map(pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI")
plot_map(pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE")
plot_map(pred_list.ZIBeta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+RE")
dev.off()