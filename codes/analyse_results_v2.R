###########################################################################
### SCRIPT TO ANALYSE RESULTS FROM FITTED BETA-REGRESSION MODELS        ###
### INCLUDES CALCULATING LOO-CV, PLOTTING MAPS, DRAWING RESPONSE CURVES ###
###########################################################################

### models fitted with scripts base.R, ZI.R, spat.R, ZI_spat.R
### this script loads fitted models and draws different maps and curves, as well as produces tables of results
library(terra)
library(loo)

# load in helper functions
source("codes/helpers.R")

# load in the training data (n = 500)
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
df_sub <- train_n500

# load in the training data (n = 1000)
#load("data/estonia_new/train/train_2020_2021_n1000.Rdata")
#df_sub <- train_n1000

# load in the training data (n = 2000)
#load("data/estonia_new/train/train_2020_2021_n2000.Rdata")
#df_sub <- train_n2000

colnames(df_sub)
colSums(df_sub[,20:38] > 0)

# remove species that are too rare (under 5 observations in this case)
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd) # see thesis for justification
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# scale the covariaes
X.scaled <- scale_covariates(X)
# add second order terms for bell-shaped response curves
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

# load in the predictive grid
predictive_grid <- vect("data/estonia_new/predictive_grid_1km_all_variables_2021_july/predictive_grid_1km_all_variables_2021_july.shp")
load("data/estonia_new/predictive_grid_1km_all_variables_2021_july.Rdata")
dim(pred_grid_1km_2021_july_df)
colnames(pred_grid_1km_2021_july_df)

# add relative light level
pred_grid_1km_2021_july_df$light_bottom <- exp(-1.7*pred_grid_1km_2021_july_df$depth/pred_grid_1km_2021_july_df$zsd)

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

# turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000


############### START ANALYZING THE RESULTS ####################

### calculate loo-CV

### MODEL 1 (BASE)
loo.beta <- c()

sp_names <- colnames(train)[20:35]
subfolder <- paste0("n_",nrow(X))

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.beta <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.beta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M1/",sp_name_modified,".rds"))
  #mod.beta_rho_min2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M1/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.beta <- c(calc_loo(mod.beta),calc_loo(mod.beta_rho))
  #loo.beta <- c(calc_loo(mod.beta),calc_loo(mod.beta_rho),calc_loo(mod.beta_rho_min2))
}

### MODEL 2
loo.ZIbeta <- c()

sp_names <- colnames(train)[20:35]
subfolder <- paste0("n_",nrow(X))

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.ZIbeta <- readRDS(paste0("models/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.ZIbeta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M2/",sp_name_modified,".rds"))
  #mod.ZIbeta_rho_min2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M2/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.ZIbeta <- c(calc_loo(mod.ZIbeta), calc_loo(mod.ZIbeta_rho))
  #loo.ZIbeta <- c(calc_loo(mod.ZIbeta), calc_loo(mod.ZIbeta_rho), calc_loo(mod.ZIbeta_rho_min2))
}

### MODEL 3
loo.beta_spat <- c()

sp_names <- colnames(train)[20:35]
subfolder <- paste0("n_",nrow(X))

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.beta_spat <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.beta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M3/",sp_name_modified,".rds"))
  #mod.beta_spat_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M3/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.beta_spat <- c(calc_loo(mod.beta_spat),calc_loo(mod.beta_spat_rho))
  #loo.beta_spat <- c(calc_loo(mod.beta_spat),calc_loo(mod.beta_spat_rho),calc_loo(mod.beta_spat_rho_min_2))
}

### MODEL 4
loo.ZIBeta_spat <- c()

sp_names <- colnames(train)[20:35]

subfolder <- paste0("n_",nrow(X.sec_ord))

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.ZIBeta_spat <- readRDS(paste0("models/",subfolder,"/M4/",sp_name_modified,".rds"))
  mod.ZIBeta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M4/",sp_name_modified,".rds"))
  #mod.ZIBeta_spat_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M4/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  loo.ZIBeta_spat <- c(calc_loo(mod.ZIBeta_spat), calc_loo(mod.ZIBeta_spat_rho))
  #loo.ZIBeta_spat <- c(calc_loo(mod.ZIBeta_spat), calc_loo(mod.ZIBeta_spat_rho), calc_loo(mod.ZIBeta_spat_rho_min_2))
}

### Combine results
loo_table <- cbind(loo.beta,loo.ZIbeta,loo.beta_spat,loo.ZIBeta_spat)
colnames(loo_table) <- c("base","ZI","RE","ZI+RE")
rownames(loo_table) <- c("common rho","rho modeled")
rownames(loo_table) <- c("rho","rho(x)>0","rho(x)>2")

### save the results
subfolder <- paste0("n_",nrow(X))
#loo_table <- cbind(loo_table, "%-zero" = 100*colMeans(train[,20:35] == 0))
save(loo_table,file = paste0("results/final_results/",subfolder,"/loo_table_rho_min_2_added.Rdata"))

### to get the latex format
library(xtable)
print(xtable(loo_table, type = "latex"))

### 2) PRODUCE MAPS (expected coverages, hotspots)

sp_names <- colnames(train)[20:35]

thinning <- 40 #take only every thinning:th posterior sample
hotspot_limit <- 0.8 #how large proportion of total coverage to preserve?

im_width <- 800
im_height <- 600

subfolder <- paste0("n_",nrow(X))

set.seed(123)
for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  ### load in models
  mod.beta <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.ZIbeta <- readRDS(paste0("models/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.beta_spat <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat <- readRDS(paste0("models/",subfolder,"/M4/",sp_name_modified,".rds"))

  ### load in models with rho modeled
  mod.beta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.ZIbeta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.beta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M4/",sp_name_modified,".rds"))

  ### make predictions over predictive grid
  pred_list.beta <- predict_beta_regression(mod.beta,
                                            pred_grid_1km_2021_july_df[,colnames(X)],
                                            X,thinning,1,rho_modeled = FALSE,C = 1000)
  pred_list.ZIbeta <- predict_ZI_beta_regression(mod.ZIbeta,
                                                 pred_grid_1km_2021_july_df[,colnames(X)],
                                                 X,thinning,1,rho_modeled = FALSE, C = 1000)
  pred_list.beta_spat <- predict_spatial_beta_regression(mod.beta_spat,
                                                         pred_grid_1km_2021_july_df[,colnames(X)],
                                                         X,
                                                         pred_grid_1km_2021_july_df[,c("x","y")],
                                                         grid_centers.df/1000, observed_grid_cells.df,thinning,1,
                                                         rho_modeled = FALSE, C = 1000)
  pred_list.ZIBeta_spat <- predict_spatial_ZI_beta_regression(mod.ZIBeta_spat,
                                                              pred_grid_1km_2021_july_df[,colnames(X)],
                                                              X,
                                                              pred_grid_1km_2021_july_df[,c("x","y")],
                                                              grid_centers.df/1000, observed_grid_cells.df,thinning,1,
                                                              rho_modeled = FALSE, C = 1000)


  ### same for model with rho modeled

  ### make predictions over predictive grid
  pred_list.beta_rho <- predict_beta_regression(mod.beta_rho,
                                                pred_grid_1km_2021_july_df[,colnames(X)],
                                                X,thinning,1,rho_modeled = TRUE,C = 1000)
  pred_list.ZIbeta_rho <- predict_ZI_beta_regression(mod.ZIbeta_rho,
                                                     pred_grid_1km_2021_july_df[,colnames(X)],
                                                     X,thinning,1,rho_modeled = TRUE, C = 1000)
  pred_list.beta_spat_rho <- predict_spatial_beta_regression(mod.beta_spat_rho,
                                                             pred_grid_1km_2021_july_df[,colnames(X)],
                                                             X,
                                                             pred_grid_1km_2021_july_df[,c("x","y")],
                                                             grid_centers.df/1000, observed_grid_cells.df,thinning,1,
                                                             rho_modeled = TRUE, C = 1000)
  pred_list.ZIBeta_spat_rho <- predict_spatial_ZI_beta_regression(mod.ZIBeta_spat_rho,
                                                                  pred_grid_1km_2021_july_df[,colnames(X)],
                                                                  X,
                                                                  pred_grid_1km_2021_july_df[,c("x","y")],
                                                                  grid_centers.df/1000, observed_grid_cells.df,thinning,1,
                                                                  rho_modeled = TRUE, C = 1000)
  
  ### Plot coverages
  png(paste0("plots/final_results/",subfolder,"/coverage_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE")
  plot_map(pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI")
  plot_map(pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE")
  plot_map(pred_list.ZIBeta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+RE")
  dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/coverage_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE+rho")
  plot_map(pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+rho")
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE+rho")
  plot_map(pred_list.ZIBeta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+RE+rho")
  dev.off()
  
  ### save the best model separately
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/coverage_maps/",sp_name_modified,"_M3.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE+rho")
  dev.off()
  
  ### save the base model separately
  png(paste0("plots/final_results/",subfolder,"/coverage_maps/",sp_name_modified,"_M1.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE")
  dev.off()
  
  ### Plot hotspots
  png(paste0("plots/final_results/",subfolder,"/hotspot_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "BASE",hotspot_limit)
  plot_map(pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI",hotspot_limit)
  plot_map(pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "RE",hotspot_limit)
  plot_map(pred_list.ZIBeta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI+RE",hotspot_limit)
  dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/hotspot_maps/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map(pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "BASE+rho",hotspot_limit)
  plot_map(pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI+rho",hotspot_limit)
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "RE+rho",hotspot_limit)
  plot_map(pred_list.ZIBeta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI+RE+rho",hotspot_limit)
  dev.off()
  
  ### save the best model separately
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/hotspot_maps/",sp_name_modified,"_M3.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "RE+rho",hotspot_limit)
  dev.off()
  
  ### save the base model separately
  png(paste0("plots/final_results/",subfolder,"/hotspot_maps/",sp_name_modified,"_M1.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map(pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "BASE",hotspot_limit)
  dev.off()
  
  ### Examine differences
  
  ##### Coverages
  # png(paste0("plots/final_results/",subfolder,"/difference_maps/coverages/4_differences/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(2,2))
  # plot_map_differences(pred_list.ZIbeta, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "ZI - base",hotspot_limit)
  # plot_map_differences(pred_list.beta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "RE - base",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "ZI+RE - ZI",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "ZI+RE - RE",hotspot_limit)
  # dev.off()
  
  # png(paste0("plots/final_results/",subfolder,"/difference_maps/coverages/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(1,1))
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "ZI+RE - base",hotspot_limit)
  # dev.off()
  # 
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/coverages/",sp_name_modified,"_M7_minus_M1.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
                       "Differences in expected coverages: (RE+rho) - BASE",hotspot_limit)
  dev.off()
  
  ##### Hotspot differences
  # png(paste0("plots/final_results/",subfolder,"/difference_maps/hotspots/4_differences/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(2,2))
  # plot_map_differences(pred_list.ZIbeta, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "ZI & base",hotspot_limit)
  # plot_map_differences(pred_list.beta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "RE & base",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) & ZI",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) & RE",hotspot_limit)
  # dev.off()
  
  # png(paste0("plots/final_results/",subfolder,"/difference_maps/hotspots/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(1,1))
  # plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) & base",hotspot_limit)
  # dev.off()
  
  # png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/hotspots/4_differences/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(2,2))
  # plot_map_differences(pred_list.ZIbeta_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "ZI & base",hotspot_limit)
  # plot_map_differences(pred_list.beta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "RE & base",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) & ZI",hotspot_limit)
  # plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) - RE",hotspot_limit)
  # dev.off()
  
  # png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/hotspots/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(1,1))
  # plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
  #                      "(ZI+RE) & base",hotspot_limit)
  # dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/hotspots/",sp_name_modified,"_M1_M7.png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map_differences(pred_list.beta_spat_rho, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "Hotspot agreement for RE+rho and base model",hotspot_limit)
  dev.off()
  
  
  ### Plot random effects for spatial models
  # with common rho
  png(paste0("plots/final_results/",subfolder,"/spatial_effects/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_random_effects(pred_list.ZIBeta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","ZI+RE (mu)")
  plot_random_effects(pred_list.ZIBeta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"pi","ZI+RE (pi)")
  plot_random_effects(pred_list.beta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","RE (mu)")
  dev.off()
  
  # rho modeled with covariates
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/spatial_effects/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_random_effects(pred_list.ZIBeta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","ZI+RE+rho (mu)")
  plot_random_effects(pred_list.ZIBeta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"pi","ZI+RE+rho (pi)")
  plot_random_effects(pred_list.beta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","RE+rho (mu)")
  dev.off()
  
}

### Calculate the Jaccard Indices between every model!
list_of_results <- list(pred_list.beta,pred_list.beta_rho,
                        pred_list.ZIbeta,pred_list.ZIbeta_rho,
                        pred_list.beta_spat,pred_list.beta_spat_rho,
                        pred_list.ZIBeta_spat,pred_list.ZIBeta_spat_rho)

jac_idx_mat <- matrix(0,nrow=8,ncol=8)

for(i in 1:length(list_of_results)) {
  pred_list_i <- list_of_results[[i]]
  for (j in ((i+1):length(list_of_results))) {
    print(paste0("i=",i,"j=",j))
    pred_list_j <- list_of_results[[j]]
    jac_idx_mat[i,j] <- plot_map_differences(pred_list_j, pred_list_i, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                                             "Hotspot agreement for RE+rho and base model",hotspot_limit,return=TRUE)
  }
}

diag(jac_idx_mat) <- 1
colnames(jac_idx_mat) <- rownames(jac_idx_mat) <- c("base","base+rho","ZI","ZI+rho","RE","RE+rho","ZI+RE","ZI+RE+rho")

save(jac_idx_mat, file = "results/final_results/n_500/jaccard_idx_table.Rdata")

### print for latex format
print(xtable(jac_idx_mat, type = "latex"))


### PRODUCE RESPONSE CURVES
sp_names <- colnames(train)[20:35]

thinning <- 20
#thinning <- 5

im_width <- 800
im_height <- 600

length_grid <- 200

subfolder <- paste0("n_",nrow(X.sec_ord))

set.seed(123)

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  ### load in models
  mod.beta <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.ZIbeta <- readRDS(paste0("models/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.beta_spat <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat <- readRDS(paste0("models/",subfolder,"/M4/",sp_name_modified,".rds"))
  
  mod_list <- list(mod.beta,mod.ZIbeta,mod.beta_spat,mod.ZIBeta_spat)
  
  ### load in rho models
  mod.beta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.ZIbeta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.beta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M4/",sp_name_modified,".rds"))
  
  mod_rho_list <- list(mod.beta_rho,mod.ZIbeta_rho,mod.beta_spat_rho,mod.ZIBeta_spat_rho)
  
  plot_and_save_responses_v2(mod_list,mod_rho_list,X,length_grid,im_width,im_height,thinning,use_median = FALSE,sp_name)
}

set.seed(123)

for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  

  ### load in rho models
  mod.beta_rho <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M1/",sp_name_modified,".rds"))
  mod.ZIbeta_rho <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M2/",sp_name_modified,".rds"))
  mod.beta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M3/",sp_name_modified,".rds"))
  mod.ZIBeta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M4/",sp_name_modified,".rds"))
  
  mod_rho_list <- list(mod.beta_rho,mod.ZIbeta_rho,mod.beta_spat_rho,mod.ZIBeta_spat_rho)
  
  plot_and_save_responses_rho_min2(mod_rho_list,X,length_grid,im_width,im_height,thinning,use_median = FALSE,sp_name)
}


### produce curves for rho
plot_and_save_rho_curves <- function(mod_list,X,grid_length = 200,im_width,im_height,thinning = 20,C=1000,use_median=FALSE,sp_name,min_rho=0) {
  # modify species name suitable for folder/file name
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  
  # create a matrix with each column capturing grid of covariate values from min to max observed values
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  subfolder <- paste0("n_",nrow(X.sec_ord))
  
  png(paste0("plots/final_results/scaled_sigmoid/rho_min_2/",subfolder,"/rho_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(mod_list))
  
  #col_names <- c("depth","nitrate","oxygen","phosphate","temperature","salinity","current","chlorophyll","light level")
  
  for(i in 1:ncol(X)) {
    print(i)
    plot(NULL, xlim = c(min(grid_matrix[,i]),max(grid_matrix[,i])),ylim=c(0,C),xlab=colnames(X)[i],ylab="rho",main="")
    #plot(NULL, xlim = c(min(grid_matrix[,i]),max(grid_matrix[,i])),ylim=c(0,C),xlab=col_names[i],ylab="rho",main="")
    mod_idx <- 1
    for(mod in mod_list) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=TRUE,1000,min_rho=min_rho)
      
      if (use_median) {
        res_rho <- apply(res$rho_sample,2,median)
      } else {
        res_rho <- colMeans(res$rho_sample)
      }
      
      #print(res_rho[1:10])
      
      lines(grid_matrix[,i],res_rho,col=cols[mod_idx])
      if (i == 1) {
        legend("topright",legend=c("base","ZI","RE","ZI+RE"),col=cols,lty=1,lwd=1)
      }
      
      mod_idx <- mod_idx + 1
    }
  }
  dev.off()
}

thinning <- 20

im_width <- 800
im_height <- 600

length_grid <- 200

subfolder <- paste0("n_",nrow(X.sec_ord))

sp_name <- sp_names[4]
sp_name_modified <- gsub(" ","_",sp_name)
sp_name_modified <- gsub("/","_",sp_name_modified)

mod.beta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M1/",sp_name_modified,".rds"))
mod.ZIbeta_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M2/",sp_name_modified,".rds"))
mod.beta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M3/",sp_name_modified,".rds"))
mod.ZIBeta_spat_rho <- readRDS(paste0("models/scaled_sigmoid/",subfolder,"/M4/",sp_name_modified,".rds"))

mod_rho_list <- list(mod.beta_rho,mod.ZIbeta_rho,mod.beta_spat_rho,mod.ZIBeta_spat_rho)

set.seed(123)
plot_and_save_rho_curves(mod_rho_list,X,length_grid,im_width,im_height,thinning,C=1000,use_median=FALSE,sp_name=sp_name,min_rho=0)


### same for rho > 2

mod.beta_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M1/",sp_name_modified,".rds"))
mod.ZIbeta_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M2/",sp_name_modified,".rds"))
mod.beta_spat_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M3/",sp_name_modified,".rds"))
mod.ZIBeta_spat_rho_min_2 <- readRDS(paste0("models/scaled_sigmoid/rho_min_2/",subfolder,"/M4/",sp_name_modified,".rds"))

mod_rho_list <- list(mod.beta_rho_min_2,mod.ZIbeta_rho_min_2,mod.beta_spat_rho_min_2,mod.ZIBeta_spat_rho_min_2)

set.seed(123)
plot_and_save_rho_curves(mod_rho_list,X,length_grid,im_width,im_height,thinning,C=1000,use_median=FALSE,sp_name=sp_name,min_rho=2)

