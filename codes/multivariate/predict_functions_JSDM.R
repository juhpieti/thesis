#######################################################################################################
### THIS SCRIPT INCLUDES PREDICTION FUNCTIONS FOR MULTIVARIATE LEFT-CENSORED BETA REGRESSION MODELS ###
### ON TOP OF BASE MODEL, ZERO-INFLATED, SPATIAL AND SPATIAL+ZERO-INFLATED ARE INCLUDED             ###
#######################################################################################################

# load in utility/helper functions
source("codes/helpers.R")

predict_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, sp_name_list, n_factors, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, b = 0.5, right_censored = FALSE, min_rho=0) {
  ### function to make predictions with multivariate left-censored beta regression
  # X.pred: prediction matrix (mxp), m locations with p covariates
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # sp_name_list: list of the names of species modeled
  # n_factors: number of latent factors
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  # min_rho: minimum possible value for rho(x) (if rho is modeled with covariates)
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
  
  ### sample latent factors for prediction locations
  # Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
  
  # initialize the list including predictions for multiple species (each as list of matrices)
  res_list <- lapply(sp_name_list, function(a) {
    list(f_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         y_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         rho_sample = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list
  
  ### loop over posterior draws
  row_idx <- 1
  for (i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
    
    ### go through J species
    for (j in 1:n_species) {
      # extract the i:th posterior sample for j:th species
      alpha_i <- post.samples$alpha[i,j]
      beta_i <- c(post.samples$beta_1[i,,j], post.samples$beta_2[i,,j])
      lambda_i <- post.samples$Lambda[i,,j]
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho_i <- post.samples$alpha_rho[i,j]
        beta_rho_i <- c(post.samples$beta_rho_1[i,,j], post.samples$beta_rho_2[i,,j])
        # rho(x) modeled with covariates
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_i + Xpred %*% beta_rho_i))
      } else {
        # every location has common rho
        rho_i <- rep(post.samples$rho[i,j], n_pred)
      }
      
      # latent f (linear predictor a + xb + zl)
      # the last zl is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)

      # mean for beta distribution using logit-link
      mu_i <- inv_logit(f_i)
      
      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,rho) with common a
      EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
      probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))
      
      ### also predict Ys
      # sample latent Vs
      V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
      
      # scale & censor
      if (!right_censored) {
        # scale to (-a,1), everything < 0 are treated as 0
        y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
      } else {
        # if right-censoring is used, scale to (-a,1+b), everything < 0 are treated as 0, everything > 1 are treated as 1
        y_i <- sapply(V_i, function(v) (max(0, min(1,(a+b+1)*v - a))))
      }
      
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i
      res_list[[j]]$EY_sam[row_idx, ] <- EYi
      res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$rho_sample[row_idx, ] <- rho_i
      res_list[[j]]$y_sam[row_idx, ] <- y_i
      
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
load("data/estonia_new/train/train_2020_2021_n100.Rdata")
X <- train_n100[,11:19]
X$light_bottom <- exp(-1.7*X$depth / X$zsd)
X <- X[,-which(colnames(X) == "zsd")] #remove secchi depth since it is not interesting for modeling in itself

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
# load in single model for Cladophora Glomerata
mod.SDM_cladophora <- readRDS("models/n_100/M1/Cladophora_glomerata.rds")
# predict
SDM_predictions <- predict_beta_regression(mod.SDM_cladophora,pred_grid_1km_2021_july_df[,colnames(X)],X,100,1,FALSE,1000,0,FALSE,0)
  
# load in JSDM model for 4 species (including Cladophora Glomerata)
mod.JSDM_4species <- readRDS("models/multivariate/n_100/9_covariates/JSDM_test.RDS")
# predict
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
JSDM_predictions <- predict_beta_regression_JSDM(mod.JSDM_4species,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,100,1,FALSE,1000,0,FALSE,0)


predict_ZI_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, sp_name_list, n_factors, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, b = 0.5, right_censored = FALSE, min_rho = 0) {
  ### function to make predictions with multivariate zero-inflated left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # sp_name_list: list of the names of species modeled
  # n_factors: number of latent factors
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # b: right-censoring constant
  # right_censored: Boolean to tell whether right-censoring is used
  # min_rho: minimum possible value for rho(x) (if rho is modeled with covariates)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) probability of suitability 
  ### 6) expectations given suitability
  ### 7) rhos
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # from df to matrix for matrix calculation
  
  # number of species and prediction locations
  n_species <- length(sp_name_list)
  n_pred <- nrow(Xpred)
  
  # extract posterior draws (all model parameters)
  post.samples <- extract(stan_fit)
  n_post_samples <- length(post.samples$lp__)
  
  # thin the posterior sample (take only every 10th)
  idx <- seq(1,n_post_samples,thinning)
  
  ### sample latent factors for prediction locations
  # Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
  # Z_pi_pred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) 
  
  # initialize the list including predictions for multiple species (each as list of matrices)
  res_list <- lapply(sp_name_list, function(a) {
    list(f_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         y_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         prob_suit_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_if_suitable_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         rho_sample = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list
  
  ### loop over posterior draws
  row_idx <- 1
  
  for(i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
    Z_pi_pred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors)
    
    ### go through J species
    for (j in 1:n_species) {
      # extract the i:th posterior sample for j:th species
      alpha_i <- post.samples$alpha[i,j]
      alpha_pi_i <- post.samples$alpha_pi[i,j]
      beta_i <- c(post.samples$beta_1[i,,j], post.samples$beta_2[i,,j])
      beta_pi_i <- c(post.samples$beta_pi_1[i,,j], post.samples$beta_pi_2[i,,j])
      lambda_i <- post.samples$Lambda[i,,j]
      lambda_pi_i <- post.samples$Lambda_pi[i,,j]
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho_i <- post.samples$alpha_rho[i,j]
        beta_rho_i <- c(post.samples$beta_rho_1[i,,j], post.samples$beta_rho_2[i,,j])
        # rho(x) modeled with covariates
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_i + Xpred %*% beta_rho_i))
      } else {
        # every location has common rho
        rho_i <- rep(post.samples$rho[i,j], n_pred)
      }
      
      # latent f (linear predictor) for mean of beta
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)

      # calculate the mean of beta using logit-link
      mu_i <- inv_logit(f_i)
      
      # latent f (linear predictor) for probability of suitability
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i + Z_pi_pred %*% lambda_pi_i)
      
      # calculate probability of suitaility using logit-link
      pi_i <- inv_logit(f_pi_i)

      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
      EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
      probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      
      ### also predict Ys
      # sample location suitabilities
      Z_i <- rbinom(nrow(Xpred),1,pi_i) #0/1 vector to tell whether the location is suitable (1 => Y from beta) or unsuitable (0 => Y=0)
      # sample latent beta variables
      V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
      
      # scale & censor
      if (!right_censored) {
        # scale to (-a,1), everything < 0 are treated as 0
        y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
      } else {
        # if right-censoring is used, scale to (-a,1+b), everything < 0 are treated as 0, everything > 1 are treated as 1
        y_i <- sapply(V_i, function(v) (max(0,min(1,(a+b+1)*v - a))))
      }
      # if the location was unsuitable, automatically y=0
      y_i <- Z_i*y_i
      
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i
      res_list[[j]]$EY_sam[row_idx, ] <- EYi
      res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$rho_sample[row_idx, ] <- rho_i
      res_list[[j]]$y_sam[row_idx, ] <- y_i
      res_list[[j]]$prob_suit_sam[row_idx, ] <- pi_i
      res_list[[j]]$EY_if_suitable_sam[row_idx, ] <- EYi_if_suitable
    }
    
    #increase the row index after every species gone through
    row_idx <- row_idx + 1 
  }
  
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### Test that it works

# load in single model for Cladophora Glomerata
mod.ZI_SDM_cladophora <- readRDS("models/n_100/M2/Cladophora_glomerata.rds")
# predict
SDM_ZI_predictions <- predict_ZI_beta_regression(mod.ZI_SDM_cladophora,pred_grid_1km_2021_july_df[,colnames(X)],X,400,1,FALSE,1000,0,FALSE,0)

# load in JSDM model for 4 species (including Cladophora Glomerata)
mod.ZI_JSDM_4species <- readRDS("models/multivariate/n_100/9_covariates/JSDM_ZI_test.RDS")
# predict
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
JSDM.ZI_predictions <- predict_ZI_beta_regression_JSDM(mod.ZI_JSDM_4species,pred_grid_1km_2021_july_df[,colnames(X)],X,sp_names,2,400,1,FALSE,1000,0,FALSE,0)


predict_spatial_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, sp_name_list, n_factors, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, min_rho = 0) {
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
         probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         rho_sample = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list
  
  # add also slot for spatially correlated latent factors
  res_list$Z_mu_sam <- array(NA, dim = c(length(idx),n_pred,n_factors)) # this is an array: matrix for each latent factor separately
  
  ### loop over posterior draws
  row_idx <- 1
  
  for(i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    Zpred <- matrix(NA,nrow=n_pred,ncol=n_factors)
    
    # sample each latent factors (spatially correlated) separately
    for (k in 1:n_factors) {
      phi_obs_ik <- post.samples$phi[i,,k] #observed (included in the training data) spatial latent factors, corrsepond to locations in S.obs
      l_ik <- post.samples$l[i,k] #length-scale of the k:th latent factor
      
      # calculate covariance block-matrices
      K_pred_obs <- exp_covariance(S.pred,S.obs,1,l_ik) # magnitude s2 = 1
      K_pred <- exp_covariance(S.pred,S.pred,1,l_ik) # magnitude s2 = 1
      K_obs <- exp_covariance(S.obs,S.obs,1,l_ik) # magnitude s2 = 1
      
      # mean and covariance for predicting phi in new locations given set of parameters and observed phi
      # justification for formulas can be read from Pietilä thesis (2025)
      phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_ik)
      phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
      
      ### sample spatial random effects over the prediction locations
      # add jitter on diagonal for computational stability
      phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
      
      # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
      L <- t(chol(phi_pred_Cov))
      phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
      
      # store k:th latent factors per location
      Zpred[,k] <- as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
      # save for output
      res_list$Z_mu_sam[row_idx,,k] <- as.vector(P %*% phi_pred_i)
    }
    
    ### go through J species
    for(j in 1:n_species) {
      # extract the i:th posterior sample for j:th species
      alpha_i <- post.samples$alpha[i,j]
      beta_i <- c(post.samples$beta_1[i,,j], post.samples$beta_2[i,,j])
      lambda_i <- post.samples$Lambda[i,,j]
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho_i <- post.samples$alpha_rho[i,j]
        beta_rho_i <- c(post.samples$beta_rho_1[i,,j], post.samples$beta_rho_2[i,,j])
        # rho(x) modeled with covariates
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_i + Xpred %*% beta_rho_i))
      } else {
        # every location has common rho
        rho_i <- rep(post.samples$rho[i,j], n_pred)
      }
      
      # latent f (linear predictor)
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)
      
      # mean of the beta distribution using logit-link
      mu_i <- inv_logit(f_i)
      
      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,rho) with common a
      EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
      probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))
      
      ### also predict Ys
      # latent beta variables
      V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
      # scale & censor
      y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
      
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i
      res_list[[j]]$EY_sam[row_idx, ] <- EYi
      res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$rho_sample[row_idx, ] <- rho_i
      res_list[[j]]$y_sam[row_idx, ] <- y_i
    }
    
    #increase the row index after every species gone through
    row_idx <- row_idx + 1 
  }
  
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### Test that it works
# load in single model for Cladophora Glomerata
mod.spat_SDM_cladophora <- readRDS("models/n_100/M3/Cladophora_glomerata.rds")
# predict
SDM_spat_predictions <- predict_spatial_beta_regression(mod.spat_SDM_cladophora,pred_grid_1km_2021_july_df[,colnames(X)],X,
                                                        pred_grid_1km_2021_july_df[,c("x","y")], grid_centers.df/1000,observed_grid_cells.df,400,1,FALSE,1000)

# load in JSDM model for 4 species (including Cladophora Glomerata)
mod.spat_JSDM_4species <- readRDS("models/multivariate/n_100/9_covariates/JSDM_spatial_test.RDS")
# predict
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")
JSDM.spat_predictions <- predict_spatial_beta_regression_JSDM(mod.spat_JSDM_4species,pred_grid_1km_2021_july_df[,colnames(X)],X,
                                                              pred_grid_1km_2021_july_df[,c("x","y")],grid_centers.df/1000,observed_grid_cells.df,
                                                              sp_names,2,400,1,FALSE,1000,0)


predict_spatial_ZI_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, sp_name_list, n_factors, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000, min_rho = 0) {
  ### function to make predictions with zero-inflated spatial left-censored beta regression
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
  ### 5) probability of suitability
  ### 6) predicted spatial random effects (for mean of beta distribution)
  ### 7) predicted spatial random effects (for probability of suitability)
  ### 8) rhos
  
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
         probzero_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         prob_suit_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         EY_if_suitable_sam = matrix(NA, nrow = length(idx), ncol = n_pred),
         rho_sample = matrix(NA, nrow = length(idx), ncol = n_pred))
  })
  names(res_list) <- sp_name_list
  
  # add also slot for spatially correlated latent factors
  res_list$Z_mu_sam <- array(NA, dim = c(length(idx),n_pred,n_factors)) # this is an array: matrix for each latent factor separately
  res_list$Z_pi_sam <- array(NA, dim = c(length(idx),n_pred,n_factors))
  
  ### loop over posterior draws
  row_idx <- 1
  
  for(i in idx) {
    ### for each posterior draw, sample latent factors for prediction locations
    Zpred <- matrix(NA,nrow=n_pred,ncol=n_factors)
    Z_pi_pred <- matrix(NA,nrow=n_pred,ncol=n_factors)
    
    # sample each latent factors (spatially correlated) separately
    for (k in 1:n_factors) {
      ### First for the Z corresponding to mean mu of Beta distribution
      phi_obs_ik <- post.samples$phi[i,,k] #observed (included in the training data) spatial latent factors, corrsepond to locations in S.obs
      l_ik <- post.samples$l[i,k] #length-scale of the k:th latent factor
      
      # calculate covariance block-matrices
      K_pred_obs <- exp_covariance(S.pred,S.obs,1,l_ik) # magnitude s2 = 1
      K_pred <- exp_covariance(S.pred,S.pred,1,l_ik) # magnitude s2 = 1
      K_obs <- exp_covariance(S.obs,S.obs,1,l_ik) # magnitude s2 = 1
      
      # mean and covariance for predicting phi in new locations given set of parameters and observed phi
      # justification for formulas can be read from Pietilä thesis (2025)
      phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_ik)
      phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
      
      ### sample spatial random effects over the prediction locations
      # add jitter on diagonal for computational stability
      phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
      
      # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
      L <- t(chol(phi_pred_Cov))
      phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
      
      # store k:th latent factors per location
      Zpred[,k] <- as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
      # save for output
      res_list$Z_mu_sam[row_idx,,k] <- as.vector(P %*% phi_pred_i)
      
      ### Exactly the same process for Z corresponding to pi (probability of suitability)
      phi_pi_obs_ik <- post.samples$phi_pi[i,,k]
      l_pi_ik <- post.samples$l_pi[i,k]
      
      # calculate covariance block-matrices
      K_pred_obs <- exp_covariance(S.pred,S.obs,1,l_pi_ik) # magnitude s2 = 1
      K_pred <- exp_covariance(S.pred,S.pred,1,l_pi_ik) # magnitude s2 = 1
      K_obs <- exp_covariance(S.obs,S.obs,1,l_pi_ik) # magnitude s2 = 1
      
      # mean and covariance for predicting phi in new locations given set of parameters and observed phi
      phi_pi_pred_m <- K_pred_obs %*% solve(K_obs,phi_pi_obs_ik)
      phi_pi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
      
      # add jitter on diagonal for computational stability
      phi_pi_pred_Cov <- phi_pi_pred_Cov + 1e-08*diag(length(phi_pi_pred_m))
      
      # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
      L <- t(chol(phi_pi_pred_Cov))
      phi_pi_pred_i <- phi_pi_pred_m + L %*% rnorm(length(phi_pi_pred_m))
      
      # store k:th latent factors per location
      Z_pi_pred[,k] <- as.vector(P %*% phi_pi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
      # save for output
      res_list$Z_pi_sam[row_idx,,k] <- as.vector(P %*% phi_pi_pred_i)
      
    }
    
    ### go through J species
    for (j in 1:n_species) {
      # extract the i:th posterior sample for j:th species
      alpha_i <- post.samples$alpha[i,j]
      alpha_pi_i <- post.samples$alpha_pi[i,j]
      beta_i <- c(post.samples$beta_1[i,,j], post.samples$beta_2[i,,j])
      beta_pi_i <- c(post.samples$beta_pi_1[i,,j], post.samples$beta_pi_2[i,,j])
      lambda_i <- post.samples$Lambda[i,,j]
      lambda_pi_i <- post.samples$Lambda_pi[i,,j]
      
      # prepare rho
      if (rho_modeled) {
        alpha_rho_i <- post.samples$alpha_rho[i,j]
        beta_rho_i <- c(post.samples$beta_rho_1[i,,j], post.samples$beta_rho_2[i,,j])
        # rho(x) modeled with covariates
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_i + Xpred %*% beta_rho_i))
      } else {
        # every location has common rho
        rho_i <- rep(post.samples$rho[i,j], n_pred)
      }
      
      # latent f (linear predictor) for mean of beta
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)
      
      # calculate the mean of beta using logit-link
      mu_i <- inv_logit(f_i)
      
      # latent f (linear predictor) for probability of suitability
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i + Z_pi_pred %*% lambda_pi_i)
      
      # calculate probability of suitability using logit-link
      pi_i <- inv_logit(f_pi_i)
      
      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
      EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
      probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      
      ### also predict Ys
      # sample suitabilities of locations
      Z_i <- rbinom(nrow(Xpred),1,pi_i) #0/1 vector to tell whether the location is suitable (1 => Y from beta) or unsuitable (0 => Y=0)
      # sample latent beta variables
      V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
      # scale & censor
      y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
      # if the location was unsuitable, automatically y = 0
      y_i <- Z_i*y_i 
      
      ### store the results for j:th species
      res_list[[j]]$f_sam[row_idx, ] <- f_i
      res_list[[j]]$EY_sam[row_idx, ] <- EYi
      res_list[[j]]$probzero_sam[row_idx, ] <- probzero_i
      res_list[[j]]$rho_sample[row_idx, ] <- rho_i
      res_list[[j]]$y_sam[row_idx, ] <- y_i
      res_list[[j]]$prob_suit_sam[row_idx, ] <- pi_i
      res_list[[j]]$EY_if_suitable_sam[row_idx, ] <- EYi_if_suitable
    }
    
    #increase the row index after every species gone through
    row_idx <- row_idx + 1 
  }
  
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### Test that it works
# load in single model for Cladophora Glomerata
mod.spat_ZI_SDM_cladophora <- readRDS("models/n_100/M4/Cladophora_glomerata.rds")
# predict
SDM_spat_ZI_predictions <- predict_spatial_ZI_beta_regression(mod.spat_ZI_SDM_cladophora,pred_grid_1km_2021_july_df[,colnames(X)],X,
                                                              pred_grid_1km_2021_july_df[,c("x","y")],grid_centers.df/1000,observed_grid_cells.df,400,1,FALSE,1000)


# load in JSDM model for 4 species (including Cladophora Glomerata)
mod.spat_ZI_JSDM_4species <- readRDS("models/multivariate/n_100/9_covariates/JSDM_spatial_ZI_test.RDS")
# predict
sp_names <- c("Cladophora glomerata","Fucus vesiculosus","Mytilus trossulus","Stuckenia pectinata")

JSDM.spat_ZI_predictions <- predict_spatial_ZI_beta_regression_JSDM(mod.spat_ZI_JSDM_4species,pred_grid_1km_2021_july_df[,colnames(X)],X,
                                                                    pred_grid_1km_2021_july_df[,c("x","y")],grid_centers.df/1000,observed_grid_cells.df,
                                                                    sp_names,2,400,1,FALSE,1000,0)
