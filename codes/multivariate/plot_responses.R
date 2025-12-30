#############################################################################################################
### SCRIPT TO DRAW RESPONSE CURVES FOR AT LEAST NON-SPATIAL VERSION OF THREE MODELS (SSDM, JSDM, 2-STAGE) ###
#############################################################################################################

# load in helpers.R for functions
load("codes/helpers.R")
load("codes/multivariate/helpers_multivariate.R")

# load in the training data (to get the covariate ranges)
load("data/estonia_new/train/train_2020_2021_all_species_n500.Rdata")
train <- train_n500_all_species

# prepare the covariate matrix
colnames(train)
X <- train[,11:19]
#X$depth_to_secchi <- X$depth / X$zsd # add secchi/depth for a variable representing seafloor light level
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

# set more comfortable names
dim(X)
colnames(X) <- c("depth","nitrate","oxygen","phosphate","temperature","salinity","current","chlorophyll","light level")

plot_responses_SDM <- function(stan_fit, X, grid_length = 200, thinning = 100, scale_curves = FALSE) {
  ### X: covariate matrix in original scale
  n_species <- length(sp_names)
  
  # prepare grid matrix (range from the observed data)
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value (this will be used to create Xpred matrices)
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  # save predictions for later use
  preds <- c()
  
  # loop over covariates
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    predX <- as.data.frame(predX)
    
    # predict
    res <- predict_beta_regression(stan_fit,predX,X,thinning,a,FALSE,1000,0,FALSE,0)
    
    # save lines
    res.EY <- res$EY_sam
    EY_line <- 100*colMeans(res.EY)
    preds <- cbind(preds, EY_line)
  }
  
  # after looping over variables, plot the results
  if(scale_curves) {
    preds <- preds / max(preds)
  }

  # plot the curves
  par(mfrow = c(3,3),
      mar = c(4,4,2,0))
  for (i in 1:ncol(X)) {
    ymin <- ifelse(scale_curves, 0,min(preds[,i]))
    ymax <- ifelse(scale_curves,1,max(preds[,i]))
    plot(grid_matrix[,i],preds[,i],type="l",ylab = "E[Y]", xlab = colnames(X)[i], main = "", ylim = c(ymin,ymax))
    }
}

grid_length <- 200
thinning <- 40
scale_curves <- FALSE

# image size for saving purposes
im_width <- 800
im_height <- 600

# subfolder based on sample size
subfolder <- paste0("n_",nrow(X))

### SSDM (non-spatial)
for (sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M1/",sp_name_modified,".rds"))
  png(paste0("plots/final_results/SSDM/",subfolder,"/M1/response_curves/resp_curves_",sp_name_modified,".png"), width = im_width, height = im_height)
  plot_responses_SDM(mod,X,grid_length,thinning,scale_curves)
  dev.off()
}

### SSDM (spatial)
for (sp_name in sp_names) {
  sp_name_modified <- gsub(" ","_",sp_name)
  sp_name_modified <- gsub("/","_",sp_name_modified)
  mod <- readRDS(paste0("models/",subfolder,"/M3/",sp_name_modified,".rds"))
  png(paste0("plots/final_results/SSDM/",subfolder,"/M3/response_curves/resp_curves_",sp_name_modified,".png"), width = im_width, height = im_height)
  plot_responses_SDM(mod,X,grid_length,thinning,scale_curves)
  dev.off()
}
  
plot_responses_JSDM <- function(stan_fit, X, sp_names, grid_length = 200, thinning = 100, scale_curves = FALSE, a=1,save_plot=FALSE,save_loc="",im_width=800,im_height=600) {
  ### X: covariate matrix in original scale
  ### sp_names: vector of species names
  n_species <- length(sp_names)
  
  # prepare grid matrix (range from the observed data)
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value (this will be used to create Xpred matrices)
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  # save predictions for all species
  preds <- list()
  
  # loop over covariates
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    predX <- as.data.frame(predX)
    
    # predict
    res <- predict_beta_regression_JSDM(stan_fit,predX,X,sp_names,2,thinning,a,FALSE,1000,0,FALSE,0,without_latent_factors = TRUE)
    
    # save lines for each species to plot later
    # loop over species
    for (j in 1:n_species) {
      res.j <- res[[sp_names[j]]]$EY_sam
      EYj_line <- 100*colMeans(res.j) #put to the scale from 0 to 100
      preds[[sp_names[j]]] <- cbind(preds[[sp_names[j]]], EYj_line)
    }
  }
  
  # after looping over variables and species, plot the results
  for (j in 1:n_species) {
    preds_j <- preds[[j]]
    
    # scale if asked
    if (scale_curves) {
      preds_j <- preds_j / max(preds_j)
    }
    
    # plot the curves
    if(save_plot) {
      sp_name <- sp_names[j]
      sp_name_modified <- gsub(" ","_",sp_name)
      sp_name_modified <- gsub("/","_",sp_name_modified)
      png(paste0(save_loc,"resp_curves_",sp_name_modified,".png"),width = im_width, height = im_height)
    }
    
    par(mfrow = c(3,3),
        mar = c(4,4,2,0))
    for (i in 1:ncol(X)) {
      ymin <- ifelse(scale_curves, 0,min(preds_j[,i]))
      ymax <- ifelse(scale_curves,1,max(preds_j[,i]))
      plot(grid_matrix[,i],preds_j[,i],type="l",ylab = "E[Y]", xlab = colnames(X)[i], main = "", ylim = c(ymin,ymax))
    }
    if (save_plot) {
      dev.off()
    }
  }
}

### JSDM (non-spatial)
mod.JSDM <- readRDS(paste0("models/multivariate/",subfolder,"/M1/JSDM_hier_priors.RDS"))
save_loc <- paste0("plots/final_results/JSDM/",subfolder,"/M1/response_curves/")
plot_responses_JSDM(mod.JSDM,X,sp_names,grid_length,thinning,scale_curves,1,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)

### JSDM (spatial)
mod.JSDM.spat <- readRDS(paste0("models/multivariate/",subfolder,"/M3/JSDM_spatial_hier_priors.RDS"))
save_loc <- paste0("plots/final_results/JSDM/",subfolder,"/M3/response_curves/")
plot_responses_JSDM(mod.JSDM.spat,X,sp_names,grid_length,thinning,scale_curves,1,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)

plot_responses_2stage <- function(stan_fit, X, sp_names, grid_length = 200, thinning = 100, scale_curves = FALSE, plot_Epi = TRUE,save_plot=FALSE,save_loc="",im_width=800,im_height=600) {
  ### X: covariate matrix in original scale
  ### sp_names: vector of species names
  ### plot_Epi: TRUE = plots expected proportions (Dirichlet part), FALSE = plots expected covers, where proportions are multiplied by expected total cover
  n_species <- length(sp_names)
  
  # prepare grid matrix (range from the observed data)
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value (this will be used to create Xpred matrices)
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  # save predictions for all species
  preds <- list()
  
  # loop over covariates
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    predX <- as.data.frame(predX)
    
    # predict
    res <- predict_DirMult_NegBin_regression(stan_fit,predX,X,sp_names,2,thinning,without_latent_factors = TRUE)
    
    # save lines for each species to plot later
    # loop over species
    for (j in 1:n_species) {
      if (plot_Epi) {
        res.j <- res[[sp_names[j]]]$Epi_sam
      } else {
        res.j <- res[[sp_names[j]]]$EY_sam
      }
      EYj_line <- colMeans(res.j)
      preds[[sp_names[j]]] <- cbind(preds[[sp_names[j]]], EYj_line)
    }
  }
  
  # after looping over variables and species, plot the results
  for (j in 1:n_species) {
    preds_j <- preds[[j]]
    
    # scale if asked
    if (scale_curves) {
      preds_j <- preds_j / max(preds_j)
    }
    
    # plot the curves
    if(save_plot) {
      sp_name <- sp_names[j]
      sp_name_modified <- gsub(" ","_",sp_name)
      sp_name_modified <- gsub("/","_",sp_name_modified)
      png(paste0(save_loc,"resp_curves_",sp_name_modified,".png"),width = im_width, height = im_height)
    }
    
    par(mfrow = c(3,3),
        mar = c(4,4,2,0))
    for (i in 1:ncol(X)) {
      ymin <- ifelse(scale_curves, 0,min(preds_j[,i]))
      ymax <- ifelse(scale_curves,1,max(preds_j[,i]))
      ylab <- ifelse(plot_Epi, "E[pi]", "E[Y]")
      plot(grid_matrix[,i],preds_j[,i],type="l",ylab = "value", xlab = colnames(X)[i], main = "", ylim = c(ymin,ymax))
    }
    
    if (save_plot) {
      dev.off()
    }
  }
}

### 2-stage (non-spatial)
mod.2stage <- readRDS(paste0("models/two_stage/",subfolder,"/M1/JSDM_NegBin_DirMult_hier_priors.RDS"))
save_loc <- paste0("plots/final_results/two_stage/",subfolder,"/M1/response_curves/expected_coverage/")
plot_responses_2stage(mod.2stage,X,sp_names,grid_length,thinning,scale_curves,plot_Epi = FALSE,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)
save_loc <- paste0("plots/final_results/two_stage/",subfolder,"/M1/response_curves/expected_proportion/")
plot_responses_2stage(mod.2stage,X,sp_names,grid_length,thinning,scale_curves,plot_Epi = TRUE,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)

### 2-stage (spatial)
mod.2stage.spat <- readRDS(paste0("models/two_stage/",subfolder,"/M3/JSDM_NegBin_DirMult_spatial_hier_priors.RDS"))
save_loc <- paste0("plots/final_results/two_stage/",subfolder,"/M3/response_curves/expected_coverage/")
plot_responses_2stage(mod.2stage.spat,X,sp_names,grid_length,thinning,scale_curves,plot_Epi = FALSE,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)
save_loc <- paste0("plots/final_results/two_stage/",subfolder,"/M3/response_curves/expected_proportion/")
plot_responses_2stage(mod.2stage.spat,X,sp_names,grid_length,thinning,scale_curves,plot_Epi = TRUE,save_plot = TRUE,save_loc = save_loc,im_width = im_width, im_height = im_height)

plot_ytot_responses_2stage <- function(stan_fit, X, grid_length = 200, thinning = 100, scale_curves = FALSE) {
  ### X: covariate matrix in original scale
  
  # prepare grid matrix (range from the observed data)
  grid_matrix <- matrix(0,nrow = grid_length, ncol = ncol(X))
  for (i in 1:ncol(X)) {
    x_grid <- seq(min(X[,i]),max(X[,i]),length=grid_length)
    grid_matrix[,i] <- x_grid
  }
  colnames(grid_matrix) <- colnames(X)
  
  # create a matrix that repeats grid_length times the mean covariate value (this will be used to create Xpred matrices)
  colmeans_matrix <- matrix(rep(colMeans(X),grid_length), nrow = grid_length, byrow = TRUE)
  colnames(colmeans_matrix) <- colnames(X)
  
  # save the predictions to re-scale later
  all_vars_preds <- c()
  
  par(mfrow = c(3,3),
      mar = c(4,4,2,0))
  
  # loop over covariates
  for (i in 1:ncol(X)) {
    predX <- colmeans_matrix # take all the variables to their mean in the training data
    predX[,i] <- grid_matrix[,i] # replace one covariate by grid
    predX <- as.data.frame(predX)
    
    # predict
    res <- predict_DirMult_NegBin_regression(stan_fit,predX,X,sp_names,2,thinning)
    res.EY <- res$EYtot_sam
    
    EY_line <- colMeans(res.EY) # posterior means are plotted
    all_vars_preds <- cbind(all_vars_preds, EY_line)
  }
  
  # if we want to scale the lines between 0 and 1:
  if (scale_curves) {
    all_vars_preds <- all_vars_preds / max(all_vars_preds) # takes values from 0 to 1
  }
    
  for (i in 1:ncol(X)) {
    ymin <- ifelse(scale_curves, 0,min(all_vars_preds[,i]))
    ymax <- ifelse(scale_curves,1,max(all_vars_preds[,i]))
    plot(grid_matrix[,i],all_vars_preds[,i],type="l",ylab = "E[y_tot]", xlab = colnames(X)[i], main = "", ylim = c(ymin,ymax))
  }
}

# draw curves for total cover
png(paste0("plots/final_results/two_stage/",subfolder,"/M1/response_curves/total_coverage.png"), width = im_width, height = im_height)
plot_ytot_responses_2stage(mod.2stage,X,grid_length,thinning,scale_curves)
dev.off()

png(paste0("plots/final_results/two_stage/",subfolder,"/M3/response_curves/total_coverage.png"), width = im_width, height = im_height)
plot_ytot_responses_2stage(mod.2stage.spat,X,grid_length,thinning,scale_curves)
dev.off()
