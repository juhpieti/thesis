### models fitted with M1-M4.R scripts
### this script loads fitted models and draws different maps and curves, as well as produces tables
library(terra)
library(loo)

# load in helper functions
source("codes/helpers.R")

# load in the training data
load("data/estonia_new/train/train_2020_2021_n500.Rdata")
df_sub <- train_n500

load("data/estonia_new/train/train_2020_2021_n1000.Rdata")
df_sub <- train_n1000

load("data/estonia_new/train/train_2020_2021_n2000.Rdata")
df_sub <- train_n2000


colnames(df_sub)
colSums(df_sub[,20:38] > 0)

# remove species that are too rare (under 5 observations)
train <- df_sub # more comfortable name
train <- train[,!(colnames(train) %in% c("Furcellaria lumbricalis loose form","Tolypella nidifica","Chara tomentosa"))]

# prepare the covariate matrix
X <- train[,11:19]
#X$depth_to_secchi <- X$depth / X$zsd # add secchi/depth for a variable representing seafloor light level
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

X.scaled <- scale_covariates(X)
### add the second order terms
X.sec_ord <- add_second_order_terms(X.scaled,colnames(X.scaled))

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
  
  # calculate leave-one-out LPD
  #loo.beta <- c(loo.beta,calc_loo(mod.beta))
  loo.beta <- c(calc_loo(mod.beta),calc_loo(mod.beta_rho))
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
  
  # calculate leave-one-out LPD
  loo.ZIbeta <- c(calc_loo(loo.ZIbeta), calc_loo(mod.ZIbeta_rho))
  #print(stan_trace(mod.ZIbeta, pars = c("alpha","rho")))
  #print(summary(mod.ZIbeta, pars = c("alpha","alpha_pi","alpha_rho"))$summary)
}

### MODEL 3
loo.beta_spat <- c()

sp_names <- colnames(train)[20:35]
subfolder <- paste0("n_",nrow(X))

for(sp_name in sp_names[1:16]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.beta_spat <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  #loo.beta_spat <- c(loo.beta_spat,calc_loo(mod.beta_spat))
  
  print(summary(mod.beta_spat, pars = c("rho","alpha","l","s2_cf"))$summary)
  print(stan_trace(mod.beta_spat, pars = c("rho","alpha","l","s2_cf")))
  
  
}

### MODEL 4
loo.ZIBeta_spat <- c()

sp_names <- colnames(train)[20:35]

subfolder <- paste0("n_",nrow(X.sec_ord))


for(sp_name in sp_names[4]) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  
  mod.ZIBeta_spat <- readRDS(paste0("models/",subfolder,"/M4/",sp_name_modified,".rds"))
  
  # calculate leave-one-out LPD
  #loo.ZIBeta_spat <- c(loo.ZIBeta_spat, calc_loo(mod.ZIBeta_spat))
  print(summary(mod.ZIBeta_spat, pars = c("alpha","alpha_pi","rho","l_mu","s2_mu","l_pi","s2_pi"))$summary)
  #print(summary(mod.ZIBeta_spat, pars = "log_lik")$summary)
  print(stan_trace(mod.ZIBeta_spat, pars = c("alpha","alpha_pi","rho","l_mu","s2_mu","l_pi","s2_pi"), inc_warmup = TRUE))
  #print(stan_trace(mod.ZIBeta_spat, pars = c("beta_pi_1"), inc_warmup = TRUE))
}

### Combine results
loo_table <- cbind(loo.beta,loo.ZIbeta,loo.beta_spat,loo.ZIBeta_spat)

colnames(loo_table) <- c("base","ZI","RE","ZI+RE")
rownames(loo_table) <- sp_names

subfolder <- paste0("n_",nrow(X))

loo_table <- cbind(loo_table, "%-zero" = 100*colMeans(train[,20:35] == 0))

save(loo_table,file = paste0("results/final_results/",subfolder,"/loo_table.Rdata"))

### to get the latex output
library(xtable)
print(xtable(loo_table, type = "latex"))



### 2) PRODUCE MAPS (expected coverages, hotspots)

sp_names <- colnames(train)[20:35]

#thinning <- 40
thinning <- 200
hotspot_limit <- 0.8

im_width <- 800
im_height <- 600

subfolder <- paste0("n_",nrow(X))

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
  plot_map(pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","BASE")
  plot_map(pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI")
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","RE")
  plot_map(pred_list.ZIBeta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover","ZI+RE")
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
           "BASE",hotspot_limit)
  plot_map(pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI",hotspot_limit)
  plot_map(pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "RE",hotspot_limit)
  plot_map(pred_list.ZIBeta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
           "ZI+RE",hotspot_limit)
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
  # png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/coverages/",sp_name_modified,".png"), width = im_width, height = im_height)
  # par(mfrow = c(1,1))
  # plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",
  #                      "ZI+RE - base",hotspot_limit)
  # dev.off()
  
  ##### Hotspots
  png(paste0("plots/final_results/",subfolder,"/difference_maps/hotspots/4_differences/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map_differences(pred_list.ZIbeta, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI - base",hotspot_limit)
  plot_map_differences(pred_list.beta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "RE - base",hotspot_limit)
  plot_map_differences(pred_list.ZIBeta_spat, pred_list.ZIbeta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - ZI",hotspot_limit)
  plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta_spat, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - RE",hotspot_limit)
  dev.off()
  
  png(paste0("plots/final_results/",subfolder,"/difference_maps/hotspots/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map_differences(pred_list.ZIBeta_spat, pred_list.beta, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - base",hotspot_limit)
  dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/hotspots/4_differences/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_map_differences(pred_list.ZIbeta_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI - base",hotspot_limit)
  plot_map_differences(pred_list.beta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "RE - base",hotspot_limit)
  plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.ZIbeta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - ZI",hotspot_limit)
  plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_spat_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - RE",hotspot_limit)
  dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/difference_maps/hotspots/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(1,1))
  plot_map_differences(pred_list.ZIBeta_spat_rho, pred_list.beta_rho, pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",
                       "ZI+RE - base",hotspot_limit)
  dev.off()
  
  
  ### Plot random effects for spatial models
  png(paste0("plots/final_results/",subfolder,"/spatial_effects/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_random_effects(pred_list.ZIBeta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","ZI+RE (mu)")
  plot_random_effects(pred_list.ZIBeta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"pi","ZI+RE (pi)")
  plot_random_effects(pred_list.beta_spat,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","RE (mu)")
  dev.off()
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/spatial_effects/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(2,2))
  plot_random_effects(pred_list.ZIBeta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","ZI+RE (mu)")
  plot_random_effects(pred_list.ZIBeta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"pi","ZI+RE (pi)")
  plot_random_effects(pred_list.beta_spat_rho,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"mu","RE (mu)")
  dev.off()

}

plot_map <- function(pred_list, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  vals <- colMeans(pred_list$EY_sam)
  vars <- apply(pred_list$EY_sam,2,var)
  test_df <- as.data.frame(cbind(locs, vals, vars))
  
  ### accumulative
  test_df_ordered <- test_df[order(-test_df$vals), ]
  test_df_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df_ordered$vals)
  idx <- which(cumsum(test_df_ordered$vals) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df_ordered[idx:nrow(test_df_ordered),"hotspot"] <- 0
  
  test_vect <- vect(test_df_ordered, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "vals")
  r.vars <- rasterize(test_vect, test_rast, field = "vars")
  r.hotspot <- rasterize(test_vect, test_rast, field = "hotspot")

  if (type == "cover") {
    plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = "E[coverage]"),
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = title_chr, col = c("red","blue"),
         plg = list(title = "hotspot", legend = c("no","yes")),
         xlab = "Easting (m)", ylab = "Northing (m)")
  }
}

plot_map_differences <- function(pred_list1,pred_list2, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  vals1 <- colMeans(pred_list1$EY_sam)
  vals2 <- colMeans(pred_list2$EY_sam)
  diff_EY <- vals1-vals2
  mean_square_diff <- mean(diff_EY^2)
  
  ### produce hotspots for model 1
  test_df1 <- as.data.frame(cbind(locs, vals1))
  
  ### accumulative
  test_df1_ordered <- test_df1[order(-test_df1$vals1), ]
  test_df1_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df1_ordered$vals1)
  idx <- which(cumsum(test_df1_ordered$vals1) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df1_ordered[idx:nrow(test_df1_ordered),"hotspot"] <- 0
  test_df1$hotspot <- test_df1_ordered$hotspot[rank(-test_df1$vals1)] #find the hotspot values
  
  ### produce hotspots for model 2
  test_df2 <- as.data.frame(cbind(locs, vals2))
  
  ### accumulative
  test_df2_ordered <- test_df2[order(-test_df2$vals2), ]
  test_df2_ordered$hotspot <- 1
  
  # find when 70% of the expected value is reached
  tot_sum <- sum(test_df2_ordered$vals2)
  idx <- which(cumsum(test_df2_ordered$vals2) > hotspot_proportion*tot_sum)[1]
  
  # set all outside to be non-hotspot
  test_df2_ordered[idx:nrow(test_df2_ordered),"hotspot"] <- 0
  test_df2$hotspot <- test_df2_ordered$hotspot[rank(-test_df2$vals2)] #find the correct hotspot values
  
  ### calculate Jaccard Index
  n_both_hotspots <- sum((test_df1$hotspot == 1) & (test_df2$hotspot == 1))
  n_mod1_hotspots <- sum(test_df1$hotspot == 1)
  n_mod2_hotspots <- sum(test_df2$hotspot == 1)
  jaccard_idx <- n_both_hotspots / (n_mod1_hotspots + n_mod2_hotspots - n_both_hotspots)
  
  shared_hotspot <- test_df1$hotspot == test_df2$hotspot
  prop_agreement <- mean(shared_hotspot)
  
  test_df <- as.data.frame(cbind(locs,diff_EY,shared_hotspot))
  
  test_vect <- vect(test_df, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "diff_EY")
  r.hotspot <- rasterize(test_vect, test_rast, field = "shared_hotspot")
  
  if (type == "cover") {
    plot(r, colNA = "lightgrey", main = paste0(title_chr," (",signif(mean_square_diff,3),")"), ,plg = list(title = "difference"),
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = paste0(title_chr," (",signif(jaccard_idx,3),")"), col = c("red","forestgreen"),
         plg = list(title = "agreement", legend = c("no","yes")),
         xlab = "Easting (m)", ylab = "Northing (m)")
  }
}

plot_random_effects <- function(pred_list, locs, pred.grid.vect, type, title_chr = "") {
  ### type: "mu" or "pi"
  
  if (type == "mu") {
    vals <- colMeans(pred_list$phi_mu_pred_sam)
    vars <- apply(pred_list$phi_mu_pred_sam,2,var)
  } else if (type == "pi") {
    vals <- colMeans(pred_list$phi_pi_pred_sam)
    vars <- apply(pred_list$phi_pi_pred_sam,2,var)
  }
  test_df <- as.data.frame(cbind(locs, vals, vars))
  
  test_vect <- vect(test_df, geom=c("x","y"), crs = "EPSG:3067")
  test_rast <- rast(ext = ext(pred.grid.vect), res = 1000, crs = "EPSG:3067")
  r <- rasterize(test_vect, test_rast, field = "vals")
  r.vars <- rasterize(test_vect, test_rast, field = "vars")

  if (type == "mu") {
    title_str <- "spatial RE"
  }

  plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = paste0("spatial RE)")),
       xlab = "Easting (m)", ylab = "Northing (m)")

}


plot_and_save_responses <- function(mod_list,X,grid_length = 200,im_width,im_height, thinning = 20) {
  
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)

  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  prob.zero_list <- list() #same for probability of zero
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod in mod_list) {
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    prob.zero_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    
    ### plots for separate curves if ZI model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      post.means <- colMeans(res.EY)
      
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        pi <- colMeans(res$prob_suit_sam)
        EY_if_suitable <- colMeans(res$EY_if_suitable_sam)
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("bottomright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    # now we want to normalize such that largest posterior men gets value 1!
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat

    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
  
    
    # draw plots
    png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i])
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  png(paste0("plots/final_results/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "relative exp. coverage", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()

  ### same for probability of zero?
  png(paste0("plots/final_results/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = "value", ylab = "prob. of zero", main = colnames(X)[i], lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
  #legend("center", legend = c("base","RE","ZI","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1.5, bty = "n")
  dev.off()
}


plot_and_save_responses_v2 <- function(mod_list,mod_rho_list,X,grid_length = 200,im_width,im_height, thinning = 20, use_median = FALSE, sp_name) {
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
  
  post.mean_list <- list() #list with n_model components, each are grid_length x p matrices
  post.mean_rho_list <- list()
  prob.zero_list <- list() #same for probability of zero
  prob.zero_rho_list <- list()
  
  mod_idx <- 1
  subfolder <- paste0("n_",nrow(X))
  for(mod_i in 1:length(mod_list)) {
    mod <- mod_list[[mod_i]]
    mod_rho <- mod_rho_list[[mod_i]]
  
    # save the posterior means
    post.mean_mat <- c() #grid_length x p matrix 
    post.mean_rho_mat <- c()
    prob.zero_mat <- c()
    prob.zero_rho_mat <- c()
    
    # save all the predictions for current model
    all_vars_preds <- list() #list with p components, each n_rep x grid_length
    all_vars_preds_rho <- list()
    
    ### plots for separate curves if ZI model
    ### first model without rho modeled
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE)
      } else {
        res <- predict_ZI_beta_regression(mod,predX,X,thinning,1,rho_modeled=FALSE)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds[[colnames(X)[i]]] <- res.EY
      
      if (use_median) {
        post.means <- apply(res.EY,2,median)
      } else {
        post.means <- colMeans(res.EY)
      }

      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        if (use_median) {
          pi <- apply(res$prob_suit_sam,2,median)
          EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
        } else {
          pi <- colMeans(res$prob_suit_sam)
          EY_if_suitable <- colMeans(res$EY_if_suitable_sam)        
        }
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), ylab = "value", xlab = colnames(X)[i], main = "",type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_mat <- cbind(post.mean_mat, post.means)
      prob.zero_mat <- cbind(prob.zero_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    
    ### then the same for log_model
    if (mod_idx %in% c(2,4)) {
      png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/M",mod_idx,"/separate_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
      par(mfrow = c(3,3),
          mar = c(5,4,2,2))
    }
    
    for (i in 1:ncol(X)) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      if (mod_idx %in% c(1,3)) {
        res <- predict_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000)
      } else {
        res <- predict_ZI_beta_regression(mod_rho,predX,X,thinning,1,rho_modeled=TRUE,1000)
      }
      res.EY <- res$EY_sam # n_rep x grid_length
      all_vars_preds_rho[[colnames(X)[i]]] <- res.EY
      
      if (use_median) {
        post.means <- apply(res.EY,2,median)
      } else {
        post.means <- colMeans(res.EY)
      }
      #plot the probability of presences and expectation given presence
      if (mod_idx %in% c(2,4)) {
        
        if (use_median) {
          pi <- apply(res$prob_suit_sam,2,median)
          EY_if_suitable <- apply(res$EY_if_suitable_sam,2,median)
        } else {
          pi <- colMeans(res$prob_suit_sam)
          EY_if_suitable <- colMeans(res$EY_if_suitable_sam)          
        }
        
        plot(grid_matrix[,i],pi,col="blue",lty=1,ylim=c(0,1), main = colnames(X)[i],type="l")
        lines(grid_matrix[,i],EY_if_suitable,col="red",lty=1)
        if (i == 1){
          legend("topright", legend = c("Pr(suitable)","E[Y|suitable]"), lty = 1, col = c("blue","red"), lwd = 1, cex = 1)
        }
        # add legend
        #plot(NULL, xlim = 0:1,ylim=0:1, ylab = "", xlab = "")
        #legend("center",legend = c("prob. suitable","E[Y|suitable]"), col = c("blue","red"), lty = 1, bty = "n", cex = 1.5)
      }
      
      # save posterior means
      post.mean_rho_mat <- cbind(post.mean_rho_mat, post.means)
      prob.zero_rho_mat <- cbind(prob.zero_rho_mat, colMeans(res$probzero_sam))
    }
    
    if (mod_idx %in% c(2,4)) {
      dev.off()
    }
    
    
    # now we want to normalize such that largest posterior men gets value 1!
    # rho not modeled
    max_post_mean <- max(post.mean_mat)
    post.mean_mat <- post.mean_mat / max_post_mean #now max value is 1
    
    # rho modeled
    max_post_mean_rho <- max(post.mean_rho_mat)
    post.mean_rho_mat <- post.mean_rho_mat / max_post_mean_rho
    
    # save scaled posterior means
    post.mean_list[[mod_idx]] <- post.mean_mat
    prob.zero_list[[mod_idx]] <- prob.zero_mat
    
    post.mean_rho_list[[mod_idx]] <- post.mean_rho_mat
    prob.zero_rho_list[[mod_idx]] <- prob.zero_rho_mat
    
    # scale also all the other samples
    all_vars_preds <- lapply(all_vars_preds, function(x) (x/max_post_mean))
    all_vars_preds_rho <- lapply(all_vars_preds_rho, function(x) (x/max_post_mean_rho))
    
    
    # draw plots
    # first for common rho model
    png(paste0("plots/final_results/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    # then for rho modeled version
    png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/M",mod_idx,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
    
    par(mfrow = c(3,3),
        mar = c(5,4,2,2))
    
    for (i in 1:ncol(X)) { # index for variable
      plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
           xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "")
      
      n_rep <- nrow(res.EY)
      for (j in 1:n_rep) { # index for sample
        lines(grid_matrix[,i],all_vars_preds_rho[[i]][j,], col = "lightgrey", lwd = 0.8)
      }
      lines(grid_matrix[,i], post.mean_rho_mat[,i], col = "black", lwd = 1)
    }
    
    dev.off()
    
    mod_idx <- mod_idx + 1
  }
  
  ### then also plot the colmeans on the same plot
  cols <- rainbow(length(post.mean_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  
  ### first the base models
  png(paste0("plots/final_results/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_list[[j]][,i], col = cols[j], lwd = 1)
      #lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lty = 2)
    }
    if (i == 1){
      #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### then the rho_modeled models
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/response_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "relative exp. coverage", main = "", lwd = 1)
    for (j in 1:length(post.mean_list)) {
      lines(grid_matrix[,i],post.mean_rho_list[[j]][,i], col = cols[j], lwd = 1)
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### same for probability of zero?
  ### first base models
  png(paste0("plots/final_results/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  #cols <- c("red","royalblue","forestgreen","yellow2")
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_list[[j]][,i], col = cols[j], lwd = 1)
      #lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lty = 2) # rho modeled version
    }
    if (i == 1){
      #legend("topright", legend = c("base","ZI","RE","ZI+RE","common rho","rho modeled"), lty = c(1,1,1,1,1,2), col = c(cols,"black","black"), lwd = 1, cex = 1)
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
  
  ### then models with rho modeled
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/prob_zero_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(prob.zero_list))
  for (i in 1:ncol(X)) {
    plot(NULL, xlim = c(min(grid_matrix[,i]), max(grid_matrix[,i])), ylim = c(0,1.1),
         xlab = colnames(X)[i], ylab = "prob. of zero", main = "", lwd = 1)
    for (j in 1:length(prob.zero_list)) {
      lines(grid_matrix[,i],prob.zero_rho_list[[j]][,i], col = cols[j], lwd=1) # rho modeled version
    }
    if (i == 1){
      legend("topright", legend = c("base","ZI","RE","ZI+RE"), lty = 1, col = cols, lwd = 1, cex = 1)
    }
  }
  dev.off()
}

### PRODUCE RESPONSE CURVES

sp_names <- colnames(train)[20:35]

thinning <- 20

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


### produce curves for rho
plot_and_save_rho_curves <- function(mod_list,X,grid_length = 200,im_width,im_height,thinning = 20,C=1000,use_median=FALSE,sp_name) {
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
  
  png(paste0("plots/final_results/scaled_sigmoid/",subfolder,"/rho_curves/",sp_name_modified,".png"), width = im_width, height = im_height)
  par(mfrow = c(3,3),
      mar = c(5,4,2,2))
  
  cols <- rainbow(length(mod_list))
  
  for(i in 1:ncol(X)) {
    print(i)
    plot(NULL, xlim = c(min(grid_matrix[,i]),max(grid_matrix[,i])),ylim=c(0,C),xlab=colnames(X)[i],ylab="rho",main="")
    mod_idx <- 1
    for(mod in mod_list) {
      predX <- colmeans_matrix # take all the variables to their mean in the training data
      predX[,i] <- grid_matrix[,i] # replace one covariate by grid
      predX <- as.data.frame(predX)
      res <- predict_beta_regression(mod,predX,X,thinning,1,rho_modeled=TRUE,1000)
      
      if (use_median) {
        res_rho <- apply(res$rho_sample,2,median)
      } else {
        res_rho <- colMeans(res$rho_sample)
      }
      
      print(res_rho[1:10])
      
      lines(grid_matrix[,i],res_rho,col=cols[mod_idx])
      if (i == 1) {
        legend("topright",legend=c("BASE","ZI","RE","ZI+RE"),col=cols,lty=1,lwd=1)
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

plot_and_save_rho_curves(mod_rho_list,X,length_grid,im_width,im_height,thinning,C=1000,use_median=FALSE,sp_name=sp_name)
