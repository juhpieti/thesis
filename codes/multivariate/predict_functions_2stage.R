###############################################################################
### THIS SCRIPT INCLUDES PREDICTION FUNCTIONS FOR 2-STAGE REGRESSION MODELS ###
### TOTAL COVERAGE IS MODELED WITH NEGATIVE-BINOMIAL                        ###
### SPECIES PROPORTIONS ARE MODELED WITH DIRICHLET                          ###
### OBSERVATION MODEL IS MULTINOMIAL                                        ###
###############################################################################

# load in utility/helper functions
source("codes/helpers.R")

library(MGLM) # for Dirichlet-Multinomial

softmax <- function(lin.preds) {
  # softmax(x) = softmax(x + C) for any C
  # substracting max(lin.preds) gives numerical stability
  Z <- lin.preds - max(lin.preds)
  return(exp(Z)/sum(exp(Z)))
}


predict_DirMult_NegBin_regression <- function(stan_fit, X.pred, X.orig, sp_name_list, n_factors, thinning = 10, without_latent_factors = FALSE) {
  ### function to make predictions with multivariate left-censored beta regression
  # X.pred: prediction matrix (mxp), m locations with p covariates
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # sp_name_list: list of the names of species modeled
  # n_factors: number of latent factors
  # thinning: use every thinning:th posterior sample in predictions
  # without_latent_factors: if TRUE, Z_i are automatically set 0, mainly due plotting purposes
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) rhos
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # df to matrix for matrix calculations
  
  # number of species and prediction locations
  n_species <- length(sp_name_list)
  n_pred <- nrow(Xpred)
  
  # extract posterior draws (all model parameters)
  post.samples <- extract(stan_fit)
  n_post_samples <- length(post.samples$lp__)
  
  # thin the posterior sample (take only every 10th)
  idx <- seq(1,n_post_samples,thinning)
  
  # initialize the list including predictions for multiple species (each as list of matrices)
  res_list <- lapply(sp_name_list, function(a) {
    list(f_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         y_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         Epi_sam = matrix(NA, nrow = length(idx), ncol = n_pred))#,
         #probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list

  # add also expected total coverages
  res_list$EYtot_sam <- matrix(NA, nrow = length(idx), ncol = n_pred) # this is an array: matrix for each latent factor separately
  res_list$Ytot_sam <- matrix(NA, nrow = length(idx), ncol = n_pred)
  
  ### loop over posterior draws
  row_idx <- 1
  for (i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    if (without_latent_factors) {
      Zpred <- matrix(0,nrow = nrow(X.pred), ncol = n_factors)
    } else {
      Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
    }

    ### Total coverage
    alpha_M_i <- post.samples$alpha_M[i]
    beta_M_i <- c(post.samples$beta_M_1[i,],post.samples$beta_M_2[i,])
    scale_M_i <- post.samples$scale_NegBin[i]
      
    fM_i <- as.vector(alpha_M_i + Xpred %*% beta_M_i)
    muM_i <- exp(fM_i)
    
    ### save the expectated total cover
    res_list$EYtot_sam[row_idx, ] <- muM_i
    
    ### Proportions
    alpha_i <- post.samples$alpha[i,] # J-vector
    beta_i <- rbind(post.samples$beta_1[i,,],post.samples$beta_2[i,,]) # PxJ matrix
    lambda_i <- post.samples$Lambda[i,,] # N_f x J matrix
    scale_i <- post.samples$scale_Dir[i] # scalar
    
    f_i <- alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i ### n_pred x J matrix
    mu_i <- t(apply(f_i,1,softmax)) ### n_pred x J matrix
    
    ### calculate expectations and probabilities of zero
    EYi <- muM_i*mu_i # expected total coverage * expected proportions
    #probzero_i <- 
    
    ### sample Ys
    M_i <- rnbinom(n_pred,size=scale_M_i,mu=muM_i) # total coverage from Negative-Binomial
    res_list$Ytot_sam[row_idx, ] <- M_i #save the prediction
    
    Y_i <- t(sapply(1:length(M_i), function(i) rdirmn(1,M_i[i],scale_i*mu_i[i,]))) #Y_i from Dirichlet-Multinomial
    
    ### store the results
    ### go through J species
    for (j in 1:n_species) {
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i[ ,j]
      res_list[[j]]$EY_sam[row_idx, ] <- EYi[ ,j]
      res_list[[j]]$Epi_sam[row_idx, ] <- mu_i[ ,j]
      #res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$y_sam[row_idx, ] <- Y_i[,j]
    }
    #increase the row index after every species gone through
    row_idx <- row_idx + 1 
  }
  
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### TEST: does it work?
## Prepare data first
# load in training data
load("data/estonia_new/train/train_2020_2021_all_species_n100.Rdata")
train <- train_n100_all_species
X <- train[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

Y <- train[,20:71]
idx_positive <- colSums(Y > 0) > 0 # take only species that have some observations
sum(idx_positive) # it's 21 of those
Y_21 <- Y[,idx_positive]
Y <- Y_21[,!colnames(Y_21) == "Ranunculus peltatus subsp_ Baudotii"] # drop 1 species for convenient J = 20

# load in the predictive grid
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
train <- train_n100
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


### Load in models and predict

# load in JSDM model for 20 species?
mod.2stage <- readRDS("models/two_stage/n_100/M1/JSDM_NegBin_DirMult.RDS")
# predict
sp_names <- colnames(Y)
JSDM_2stage_predictions <- predict_DirMult_NegBin_regression(mod.2stage,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,100)

predict_spatial_DirMult_NegBin_regression <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, sp_name_list, n_factors, thinning = 10) {
  ### function to make predictions with spatial left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # pred.locations: locations for each prediction cell (in meters)
  # S.pred: locations of coarse spatial grid cells (in kilometers)
  # S.obs: locations of observed coarse spatial grid cells (in kilometers, WATCH THAT THIS IS THE SAME AS WHEN FITTING A MODEL)
  # sp_name_list: list of the names of species modeled
  # n_factors: number of latent factors
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # min_rho: minimum possible value for rho(x) (if rho is modeled with covariates)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) predicted spatial random effects
  ### 6) rhos
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # df to matrix for matrix calculations
  
  # prepare P matrix (mxS) that tells in which spatial random effect cell each prediction point belongs to
  P <- prepare_P_matrix(pred.locations,S.pred*1000) # turn grid locations to meters
  
  # number of species and prediction locations
  n_species <- length(sp_name_list)
  n_pred <- nrow(Xpred)
  
  # extract posterior draws (all model parameters)
  post.samples <- extract(stan_fit)
  n_post_samples <- length(post.samples$lp__)
  
  # thin the posterior sample (take only every 10th)
  idx <- seq(1,n_post_samples,thinning)
  
  # initialize the list including predictions for multiple species (each as list of matrices)
  res_list <- lapply(sp_name_list, function(a) {
    list(f_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         y_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         Epi_sam = matrix(NA, nrow = length(idx), ncol = n_pred))#,
    #probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list
  
  # add also expected total coverages
  res_list$EYtot_sam <- matrix(NA, nrow = length(idx), ncol = n_pred) # this is an array: matrix for each latent factor separately
  res_list$Ytot_sam <- matrix(NA, nrow = length(idx), ncol = n_pred)
  
  # add also slot for spatially correlated latent factors
  res_list$Z_mu_sam <- array(NA, dim = c(length(idx),n_pred,n_factors)) # this is an array: matrix for each latent factor separately
  
  # and spatially correlated random effects for Neg-Bin part of the model (total coverage)
  res_list$phi_M_sam <- matrix(NA, nrow = length(idx), ncol = n_pred)
  
  ### loop over posterior draws
  row_idx <- 1
  
  for(i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    #Zpred <- matrix(NA,nrow=n_pred,ncol=n_factors)
    # NOTE, SPATIALLY CORRELATED PART WILL BE ADDED LATER ON
    Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,0.5),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
    
    # sample each latent factors (spatially correlated) separately
    for (k in 1:n_factors) {
      phi_obs_ik <- post.samples$phi[i,,k] #observed (included in the training data) spatial latent factors, corresponds to locations in S.obs
      l_ik <- post.samples$l[i,k] #length-scale of the k:th latent factor
      
      # calculate covariance block-matrices
      K_pred_obs <- exp_covariance(S.pred,S.obs,0.5,l_ik) # magnitude s2 = 0.5
      K_pred <- exp_covariance(S.pred,S.pred,0.5,l_ik) # magnitude s2 = 0.5
      K_obs <- exp_covariance(S.obs,S.obs,0.5,l_ik) # magnitude s2 = 0.5
      
      # mean and covariance for predicting phi in new locations given set of parameters and observed phi
      # justification for formulas can be read from Pietilä thesis (2025)
      phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_ik)
      phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
      
      ### sample spatial random effects over the prediction locations
      # add jitter on diagonal for computational stability
      phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))

      # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
      L <- t(chol(phi_pred_Cov))
      #print(paste0(i,k,"first cholesky"))
      phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
      
      # store k:th latent factors per location
      Zpred[,k] <- Zpred[,k] + as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
      # save for output
      res_list$Z_mu_sam[row_idx,,k] <- as.vector(P %*% phi_pred_i)
    }
    
    ### sample spatial random effects for Neg-Bin model
    phi_M_obs_i <- post.samples$phi_M[i, ] #observed random effects
    l_M_i <- post.samples$l_M[i] #length-scale
    s2_M_i <- post.samples$s2_M[i] #magnitude
    
    # calculate covariance block-matrices
    KM_pred_obs <- exp_covariance(S.pred,S.obs,s2_M_i,l_M_i)
    KM_pred <- exp_covariance(S.pred,S.pred,s2_M_i,l_M_i) 
    KM_obs <- exp_covariance(S.obs,S.obs,s2_M_i,l_M_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_M_pred_m <- KM_pred_obs %*% solve(KM_obs,phi_M_obs_i)
    phi_M_pred_Cov <- KM_pred - KM_pred_obs %*% solve(KM_obs) %*% t(KM_pred_obs)
    
    ### sample spatial random effects over the prediction locations
    # add jitter on diagonal for computational stability
    phi_M_pred_Cov <- phi_M_pred_Cov + 1e-08*diag(length(phi_M_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    LM <- t(chol(phi_M_pred_Cov))
    phi_M_pred_i <- phi_M_pred_m + LM %*% rnorm(length(phi_M_pred_m))
    
    res_list$phi_M_sam[row_idx, ] <- as.vector(P %*% phi_M_pred_i)
    
    ### Total coverage
    alpha_M_i <- post.samples$alpha_M[i]
    beta_M_i <- c(post.samples$beta_M_1[i,],post.samples$beta_M_2[i,])
    scale_M_i <- post.samples$scale_NegBin[i]
    
    fM_i <- as.vector(alpha_M_i + Xpred %*% beta_M_i + P %*% phi_M_pred_i) # fixed effects + random effects
    muM_i <- exp(fM_i)
    
    ### save the expectated total cover
    res_list$EYtot_sam[row_idx, ] <- muM_i
    
    ### Proportions
    alpha_i <- post.samples$alpha[i,] # J-vector
    beta_i <- rbind(post.samples$beta_1[i,,],post.samples$beta_2[i,,]) # PxJ matrix
    lambda_i <- post.samples$Lambda[i,,] # N_f x J matrix
    scale_i <- post.samples$scale_Dir[i] # scalar
    
    f_i <- alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i ### n_pred x J matrix
    mu_i <- t(apply(f_i,1,softmax)) ### n_pred x J matrix
    
    ### calculate expectations and probabilities of zero
    EYi <- muM_i*mu_i # expected total coverage * expected proportions
    #probzero_i <- 
    
    ### sample Ys
    M_i <- rnbinom(n_pred,size=scale_M_i,mu=muM_i) # total coverage from Negative-Binomial
    res_list$Ytot_sam[row_idx, ] <- M_i #save the prediction
    
    Y_i <- t(sapply(1:length(M_i), function(i) rdirmn(1,M_i[i],scale_i*mu_i[i,]))) #Y_i from Dirichlet-Multinomial
    
    ### store the results
    ### go through J species
    for (j in 1:n_species) {
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i[ ,j]
      res_list[[j]]$EY_sam[row_idx, ] <- EYi[ ,j]
      res_list[[j]]$Epi_sam[row_idx, ] <- mu_i[ ,j]
      #res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$y_sam[row_idx, ] <- Y_i[,j]
    }
    
    #increase the row index after every species gone through
    row_idx <- row_idx + 1 
  }
  
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### Test that it works
# load in JSDM model for 20 species?
mod.2stage.spat <- readRDS("models/two_stage/n_100/M3/JSDM_NegBin_DirMult_spatial.RDS")
# predict
sp_names <- colnames(Y)
JSDM_2stage_spatial_predictions <- predict_spatial_DirMult_NegBin_regression(mod.2stage.spat,pred_grid_1km_2021_july_df[,colnames(X)],X,
                                                                             pred_grid_1km_2021_july_df[,c("x","y")],grid_centers.df/1000,observed_grid_cells.df,
                                                                             sp_names,2,100)



