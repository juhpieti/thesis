library(xtable) # to output latex tables
library(terra) # spatial prediction functions use this package

# load in helper functions
load("codes/helpers.R")

### bring in the training data
load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
dim(train_n500_all_species)
Xtrain <- train_n500_all_species[,11:19]
Xtrain$light_bottom <- exp(-1.7*Xtrain$depth / Xtrain$zsd)
Xtrain <- Xtrain[,-which(colnames(Xtrain) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

Ytrain <- train_n500_all_species[,20:71]
idx_positive <- colSums(Ytrain > 0) > 2
sum(idx_positive)
Ytrain <- Ytrain[,idx_positive]
Ytrain <- Ytrain[,!colnames(Ytrain) == "Ranunculus peltatus subsp_ Baudotii"] # drop 1 species for convenient J = 20
dim(Ytrain)

### bring in the test data
load("data/estonia_new/test/test_2020_2021_all_species_n500.Rdata")
dim(test_n500_all_species)

Xtest <- test_n500_all_species[,11:19]
Xtest$light_bottom <- exp(-1.7*Xtest$depth / Xtest$zsd)
Xtest <- Xtest[,-which(colnames(Xtest) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

sp_names <- colnames(Ytrain)
Ytest <- test_n500_all_species[,sp_names]
dim(Ytest)

### prepare the data for spatial models
### load in the coarse spatial random effect grid
spatial_grid <- vect("data/estonia_new/spatial_random_effect_grid_20km/spatial_random_effect_grid_20km.shp")

# prepare a matrix of grid centers and their locations (TM35FIN coordinates)
grid_centers <- centroids(spatial_grid)
grid_centers.df <- as.data.frame(grid_centers, geom = "XY")
grid_centers.df <- grid_centers.df[,c("x","y")]
dim(grid_centers.df) # there are m = 191 grid cells

# turn the training data into terra-vector
train <- train_n500_all_species
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

# turn the coordinates in km instead of meters
observed_grid_cells.df <- observed_grid_cells.df/1000

### make predictions for test set
thinning <- 40 #do not use all of the posterior sample

# STACKED SDMS (non-spatial)
n_species <- length(sp_names)
test_pred_SSDM <- list()
subfolder <- paste0("n_",nrow(Xtrain)) # find the correct folder for models

for (sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  pred <- predict_beta_regression(mod,Xtest,Xtrain,thinning,1,FALSE,1000,0,FALSE,0)
  test_pred_SSDM[[sp_name]] <- pred
}

# STACKED SDMS (spatial)
test_pred_SSDM.spat <- list()
for (sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  pred <- predict_spatial_beta_regression(mod,Xtest,Xtrain,test_n500_all_species[,c("x","y")],
                                          grid_centers.df/1000,observed_grid_cells.df,thinning,1,FALSE,1000)
  test_pred_SSDM.spat[[sp_name]] <- pred
}

# JSDM (non-spatial)
mod.JSDM <- readRDS("models/multivariate/n_500/M1/JSDM_hier_priors.RDS")
set.seed(123)
test_pred_JSDM <- predict_beta_regression_JSDM(mod.JSDM, Xtest, Xtrain, sp_names, 2, thinning, 1, FALSE, 1000, 0, FALSE, 0)

# JSDM (spatial)
mod.JSDM.spat <- readRDS("models/multivariate/n_500/M3/JSDM_spatial_hier_priors.RDS")
set.seed(123)
test_pred_JSDM.spat <- predict_spatial_beta_regression_JSDM(mod.JSDM.spat, Xtest, Xtrain, test_n500_all_species[,c("x","y")],
                                                            grid_centers.df/1000,observed_grid_cells.df,sp_names,2,thinning,1,FALSE,1000,0)


# 2-STAGE MODEL (non-spatial)
mod.2stage <- readRDS("models/two_stage/n_500/M1/JSDM_NegBin_DirMult_hier_priors.RDS")
set.seed(123)
test_pred_2stage <- predict_DirMult_NegBin_regression(mod.2stage, Xtest, Xtrain, sp_names, 2, thinning)

# 2-STAGE MODEL (spatial)
mod.2stage.spat <- readRDS("models/two_stage/n_500/M3/JSDM_NegBin_DirMult_spatial_hier_priors.RDS")
set.seed(123)
test_pred_2stage.spat <- predict_spatial_DirMult_NegBin_regression(mod.2stage.spat,Xtest,Xtrain,test_n500_all_species[,c("x","y")],
                                                                   grid_centers.df/1000,observed_grid_cells.df,sp_names,2,thinning)

### write functions
scale_predictions <- function(pred_list, sp_names, scaling_coeff = 100) {
  ### scale predictions with 100 or 0.01 (or anything else)
  pred_list_new <- pred_list
  for (sp_name in sp_names) {
    pred_list_new[[sp_name]]$y_sam <- scaling_coeff*pred_list[[sp_name]]$y_sam
    pred_list_new[[sp_name]]$EY_sam <- scaling_coeff*pred_list[[sp_name]]$EY_sam
  }
  return(pred_list_new)
}

### THE CODE STARTING FROM HERE IS UGLY AND WRITTEN IN RUSH! REWRITTEN IF TIME LEFT

compare_predictions <- function(pred_list, Y_obs, use_expectation = FALSE, use_median = FALSE) {
  ### pred_list: posterior predictions for test set
  ### Y_obs: observed Y values (test set)
  ### use_expectation: instead of sampled Y, use E[Y] to compare to observed Y
  ### use_median: instead of mean, use median to summarize the posterior distribution
  
  # calculate observed proportions
  sum_obs <- rowSums(Y_obs)
  P_obs <- Y_obs / (sum_obs + 1e-12)
  
  # initialize list for results
  res_list <- list()
  
  # number of species
  n_species <- ncol(Y_obs)
  
  # number of potserior draws
  n_rep <- nrow(pred_list[[1]][[1]])
  
  # initialize lists for results
  dists_sp_y <- matrix(0, nrow = n_rep, ncol = n_species)
  dists_sp_p <- matrix(0, nrow = n_rep, ncol = n_species)
  
  dists_community_Y <- c()
  dists_community_P <- c()
  
  dists_tot <- c()
  errors_tot <- c()
  
  for(i in 1:n_rep) {
    ### gather replicated data Y by looping over species
    Y_rep <- c()
    for (j in 1:n_species) {
      if (use_expectation) {
        Y_rep <- cbind(Y_rep, pred_list[[j]][["EY_sam"]][i, ])
      } else {
        Y_rep <- cbind(Y_rep, pred_list[[j]][["y_sam"]][i, ])
      }
    }

    ### calculate proportion matrix for replicated observations
    sum_rep <- rowSums(Y_rep) # total coverage
    P_rep <- Y_rep / (sum_rep + 1e-12) # proportions
    
    ### Calculate specieswise RMSE, average is taken over sampling locations
    for (j in 1:n_species) {
      dists_sp_y[i,j] <- sqrt(mean((Y_rep[,j]-Y_obs[,j])^2))
      dists_sp_p[i,j] <- sqrt(mean((P_rep[,j]-P_obs[,j])^2))
    }
    
    ### Calculate also distances ||y_i - ỹ_i|| and ||p_i - ~p_i||
    eucl_dists_Y <- sqrt(rowSums((Y_obs - Y_rep)^2))
    eucl_dists_P <- sqrt(rowSums((P_obs - P_rep)^2))
    
    ### Save rooted average over locations
    dists_community_Y <- c(dists_community_Y, sqrt(mean(eucl_dists_Y)))
    dists_community_P <- c(dists_community_P, sqrt(mean(eucl_dists_P)))
    
    ### Calculate RMSE for total coverage
    dists_tot <- c(dists_tot, sqrt(mean((sum_rep - sum_obs)^2)))
    
    ### Calculate mean difference between predicted and observed total coverage
    errors_tot <- c(errors_tot, mean(sum_rep - sum_obs))

  # if (use_median) {
  #   return(list(sp_dist_Y = apply(mean_sp_dists_y,2,median),
  #               sp_dist_P = apply(mean_sp_dists_p,2,median),
  #               eucl_dist_Y = median(mean_eucl_dists_community_Y),
  #               eucl_dist_P = median(mean_eucl_dists_community_Y),
  #               squared_error_tot = median(mean_squared_errors_tot),
  #               error_tot = median(mean_errors_tot)))
  # } else {
  #   return(list(sp_dist_Y = colMeans(mean_sp_dists_y),
  #               sp_dist_P = colMeans(mean_sp_dists_p),
  #               eucl_dist_Y = mean(mean_eucl_dists_community_Y),
  #               eucl_dist_P = mean(mean_eucl_dists_community_P),
  #               squared_error_tot = mean(mean_squared_errors_tot),
  #               error_tot = mean(mean_errors_tot)))
  }
  return(list(sp_dist_Y = dists_sp_y,
              sp_dist_P = dists_sp_p,
              comm_dist_Y = dists_community_Y,
              comm_dist_P = dists_community_P,
              dist_tot = dists_tot,
              error_tot = errors_tot))
}

### run comparisons
# initialize the results table
res_table <- c()
use_median <- TRUE #use posterior medians instead of means
use_expectation <- FALSE #use E[Y] instead of sampled Y as y_hat

### SSDM (non-spatial)
test_pred_SSDM <- scale_predictions(test_pred_SSDM, sp_names, 100)
comp_res_SSDM <- compare_predictions(test_pred_SSDM, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_SSDM$sp_dist_Y,1,mean)),
                                median(comp_res_SSDM$comm_dist_Y),
                                median(comp_res_SSDM$dist_tot),
                                median(comp_res_SSDM$error_tot,2,median)))

### SSDM (spatial)
test_pred_SSDM.spat <- scale_predictions(test_pred_SSDM.spat, sp_names, 100)
comp_res_SSDM.spat <- compare_predictions(test_pred_SSDM.spat, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_SSDM.spat$sp_dist_Y,1,mean)),
                                median(comp_res_SSDM.spat$comm_dist_Y),
                                median(comp_res_SSDM.spat$dist_tot),
                                median(comp_res_SSDM.spat$error_tot)))


### JSDM (non-spatial)
test_pred_JSDM <- scale_predictions(test_pred_JSDM, sp_names, 100) # scale to same scale with 2stage model
comp_res_JSDM <- compare_predictions(test_pred_JSDM, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_JSDM$sp_dist_Y,1,mean)),
                                median(comp_res_JSDM$comm_dist_Y),
                                median(comp_res_JSDM$dist_tot),
                                median(comp_res_JSDM$error_tot)))

### JSDM (spatial)
test_pred_JSDM.spat <- scale_predictions(test_pred_JSDM.spat, sp_names, 100)
comp_res_JSDM.spat <- compare_predictions(test_pred_JSDM.spat, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_JSDM.spat$sp_dist_Y,1,mean)),
                                median(comp_res_JSDM.spat$comm_dist_Y),
                                median(comp_res_JSDM.spat$dist_tot),
                                median(comp_res_JSDM.spat$error_tot)))


### 2-stage (non-spatial)
comp_res_2stage <- compare_predictions(test_pred_2stage, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_2stage$sp_dist_Y,1,mean)),
                                median(comp_res_2stage$comm_dist_Y),
                                median(comp_res_2stage$dist_tot),
                                median(comp_res_2stage$error_tot)))


### 2-stage (spatial)
comp_res_2stage.spat <- compare_predictions(test_pred_2stage.spat, Ytest, use_expectation, use_median)

res_table <- cbind(res_table, c(median(apply(comp_res_2stage.spat$sp_dist_Y,1,mean)),
                                median(comp_res_2stage.spat$comm_dist_Y),
                                median(comp_res_2stage.spat$dist_tot),
                                median(comp_res_2stage.spat$error_tot)))

# set the row and column names
rownames(res_table) <- c("y_ij","y_i","y_tot (MSE)", "y_tot (diff)")
colnames(res_table) <- c("SSDM","SSDM+RE","JSDM","JSDM+RE","2stage","2stage+RE")
res_table

# output the table as latex table
print(xtable(res_table, type = "latex"))

# plot posterior distributions of the predictive metrics
mod_names <- c("SSDM","SSDM+RE","JSDM","JSDM+RE","2stage","2stage+RE")
pred_metrics <- c("y_ij","p_ij","Y_i","P_i","Ytot (RMSE)", "Ytot (DIFF)")

im_width <- 800
im_height <- 600

png(paste0("plots/final_results/predictive_model_comparison_Y.png"), width = im_width, height = im_height)

par(mfrow = c(2,3))
list_of_preds <- list(comp_res_SSDM,
                      comp_res_SSDM.spat,
                      comp_res_JSDM,
                      comp_res_JSDM.spat,
                      comp_res_2stage,
                      comp_res_2stage.spat)

for(i in 1:length(pred_metrics)) {
  data_mat <- c()
  for (m in 1:length(list_of_preds)) {
    if (i %in% 1:2) {
      data_mat <- cbind(data_mat, rowMeans(list_of_preds[[m]][[i]])) #m:th model, i:th metrics
    } else {
      data_mat <- cbind(data_mat, list_of_preds[[m]][[i]]) #m:th model, i:th metrics
    }
  }
  boxplot(data_mat,names = mod_names, main = pred_metrics[i], xlab = "model", ylab = "value", cex.axis = 0.8)
}

dev.off()

# make table for species-level rank of the models
res_species <- cbind(apply(comp_res_SSDM$sp_dist_Y,2,median),
                     apply(comp_res_SSDM.spat$sp_dist_Y,2,median),
                     apply(comp_res_JSDM$sp_dist_Y,2,median),
                     apply(comp_res_JSDM.spat$sp_dist_Y,2,median),
                     apply(comp_res_2stage$sp_dist_Y,2,median),
                     apply(comp_res_2stage.spat$sp_dist_Y,2,median))

rownames(res_species) <- sp_names
colnames(res_species) <- c("SSDM","SSDM+RE","JSDM","JSDM+RE","2stage","2stage+RE")

res_species <- cbind(res_species, prevalence = 100*colMeans(Ytest > 0))
res_species

# latex table output
print(xtable(res_species, type = "latex"))

