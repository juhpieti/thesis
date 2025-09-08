###########################################################################
### SCRIPT TO ANALYSE RESULTS FROM FITTED JSDM BETA-REGRESSION MODELS   ###
### INCLUDES CALCULATING LOO-CV, PLOTTING MAPS, DRAWING RESPONSE CURVES ###
###########################################################################

### models fitted with separate R-scripts (joint_beta_regression.R, joint_zi_beta_regression.R, joint_beta_regression_spatial.R, joint_zi_beta_regression_spatial.R)
### this script loads fitted models and draws different maps and curves, as well as produces tables
library(terra)
library(loo)

# load in helper functions
source("codes/helpers.R")

# load in the training data
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
train <- train_n100

# load("data/estonia_new/train/train_2020_2021_n500.Rdata")
# train <- train_n500

colnames(train)
colSums(train[,20:38] > 0)

# remove species that are too rare (under 5 observations)
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

### prepare functions

plot_map <- function(pred_mat, locs, pred.grid.vect, type, title_chr = "", hotspot_proportion = 0.7) {
  ### type: "cover", "hotspot", 
  
  #vals <- colMeans(pred_list$EY_sam)
  vals <- colMeans(pred_mat)
  #vars <- apply(pred_list$EY_sam,2,var)
  vars <- apply(pred_mat,2,var)
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
    # plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = "E[coverage]"),
    #      xlab = "Easting (m)", ylab = "Northing (m)")
    plot(r, colNA = "lightgrey", main = title_chr,
         xlab = "Easting (m)", ylab = "Northing (m)")
  } else if (type == "hotspot") {
    plot(r.hotspot, colNA = "lightgrey", main = title_chr, col = c("red","blue"),
         #plg = list(title = "hotspot", legend = c("no","yes")),
         legend = FALSE,
         xlab = "Easting (m)", ylab = "Northing (m)")
    #legend(583000,6362000, legend = c("yes","no"), title = "hotspot", col = c("blue","red"), pch = 15,cex = 1)
    legend("bottomright", legend = c("yes","no"), title = "hotspot", col = c("blue","red"), pch = 15,cex = 1, inset = c(0.001,0.04))
  }
}

plot_map_JSDM <- function(pred_lists, locs, pred.grid.vect, type, hotspot_proportion, sp_names, summary_maps = FALSE) {
  ###
  # pred_lists: named list (by species name), each component is a list of predictions as a matrix (n_posterior_samples x n_prediction locations)
  
  n_species <- length(sp_names)
  
  if (!summary_maps) {
    # prepare a grid for plotting, based on the number of species
    n_rows <- floor(sqrt(n_species)) #prefer square type grid, e.g. 2x2, 3x3 if possible
    n_cols <- ceiling(n_species/n_rows) #add so many columns that every species fit
    par(mfrow = c(n_rows,n_cols))
    
    # plot each species separately
    for(sp_name in sp_names) {
      pred_list <- pred_lists[[sp_name]]
      plot_map(pred_list$EY_sam,locs,pred.grid.vect,type,sp_name,hotspot_proportion)
    }
  } else {
    ### 4 different plots
    ## 1) total percent cover
    ## 2) species richness
    ## 3) Simpson index
    ## 4) Shannon index
    
    # initialize matrices
    nr <- nrow(pred_lists[[1]]$EY_sam)
    nc <- ncol(pred_lists[[1]]$EY_sam)
    
    sum_mat <- sp_rich_mat <- simpson_mat <- shannon_mat <- matrix(0,nrow=nr,ncol=nc)
    
    # first calculate total percent cover and no. species present
    for (j in 1:n_species) {
      sum_mat <- sum_mat + pred_lists[[j]]$EY_sam
      sp_rich_mat <- sp_rich_mat + (pred_lists[[j]]$y_sam > 0) 
    }
    
    # for diversity indeces, calculate proportional covers p_ij, such that sum_j(p_ij) = 1 for each location i
    for (j in 1:n_species) {
      # relative cover for species j
      p_j <- pred_lists[[j]]$EY_sam / sum_mat
      # add j:th species contribution to simpson idx
      simpson_mat <- simpson_mat + p_j^2
      # add j:th species contribution to shannon idx
      shannon_mat <- shannon_mat + p_j*log(p_j + 1e-12)
    }
    
    # simpson index is defined as 1 - sum_j(p_j^2)
    simpson_mat <- 1 - simpson_mat
    
    # shannon index is defined as -1*sum_j(p_j*ln(p_j))
    shannon_mat <- -1*shannon_mat
    
    par(mfrow = c(2,2))
    plot_map(sum_mat,locs,pred.grid.vect,type,"total percent cover",hotspot_proportion)
    plot_map(sp_rich_mat,locs,pred.grid.vect,type,"species richness",hotspot_proportion)
    plot_map(simpson_mat,locs,pred.grid.vect,type,"Simpson index",hotspot_proportion)
    plot_map(shannon_mat,locs,pred.grid.vect,type,"Shannon index",hotspot_proportion)
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
  
  plot(r, colNA = "lightgrey", main = title_chr, ,plg = list(title = paste0("spatial RE")),
       xlab = "Easting (m)", ylab = "Northing (m)")
}

plot_random_effects_JSDM <- function(pred_list,locs,pred.grid.vect,type,title_chr="") {
  
}



### DRAW MAPS ###
sp_names <- names(JSDM_predictions)

##### LEFT-CENSORED BETA REGRESSION (M1)

### 1) Stacked SDM
subfolder <- paste0("n_",nrow(X)) # find the correct folder for models

# gather the predictions as lists
stacked_SDM_predictions <- list()
for(sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  pred <- predict_beta_regression(mod,pred_grid_1km_2021_july_df[,colnames(X)],X,100,1,FALSE,1000,0,FALSE,0)
  stacked_SDM_predictions[[sp_name]] <- pred
}

# plot with JSDM function
plot_map_JSDM(stacked_SDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)
plot_map_JSDM(stacked_SDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,TRUE)

# plot also hotspots
plot_map_JSDM(stacked_SDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.7,sp_names,FALSE)
plot_map_JSDM(stacked_SDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.7,sp_names,TRUE)


### 2) JSDM

# load model

# predict 

# plot
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,TRUE)

# plot also hotspots
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.7,sp_names,FALSE)
plot_map_JSDM(JSDM_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"hotspot",0.7,sp_names,TRUE)


##### 2) ZERO-INFLATED LEFT-CENSORED BETA REGRESSION
plot_map_JSDM(JSDM.ZI_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)

##### 3) SPATIAL LEFT-CENSORED BETA REGRESSION
plot_map_JSDM(JSDM.spat_predictions,pred_grid_1km_2021_july_df[,c("x","y")],predictive_grid,"cover",0.7,sp_names,FALSE)

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