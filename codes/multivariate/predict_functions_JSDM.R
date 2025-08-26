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
  
  ### sample latent factors for prediction locations
  Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
  
  # initialize the list including predictions for multiple species
  res_list <- list()
  
  ### loop over species
  n_species <- length(sp_name_list)
  
  for(j in 1:n_species) {
    # load the posterior samples (take only j:th species)
    alpha_sam <- as.matrix(stan_fit, pars = paste0("alpha[",j,"]"))
    beta_sam <- as.matrix(stan_fit, pars = c(paste0("beta_1[",1:ncol(X.pred),",",j,"]"),paste0("beta_2[",1:ncol(X.pred),",",j,"]")))
    lambda_sam <- as.matrix(stan_fit, pars = paste0("Lambda[",1:n_factors,",",j,"]"))
    
    # prepare rho
    if (rho_modeled) {
      # coefficient related to rho(x)
      alpha_rho_sam <- as.matrix(stan_fit, pars = paste0("alpha_rho[",j,"]"))
      beta_rho_sam <- as.matrix(stan_fit, pars = c(paste0("beta_rho_1[",1:ncol(X.pred),",",j,"]"),paste0("beta_rho_2[",1:ncol(X.pred),",",j,"]")))
    } else {
      # posterior sample of common rho
      rho_sam <- as.matrix(stan_fit, pars = paste0("rho[",j,"]")) #j:th species
    }

    # initialize matrices for output
    f_sam <- c()
    y_sam <- c()
    EY_sam <- c()
    probzero_sam <- c()
    rho_sample <- c()
    
    # loop over posterior samples (take every thinning:th sample)
    for (i in seq(1,nrow(beta_sam),thinning)) {
      # corresponding coefficients
      alpha_i <- alpha_sam[i,]
      beta_i <- beta_sam[i,]
      lambda_i <- lambda_sam[i,]
      
      # prepare rho
      if (rho_modeled){
        # rho(x) modeled with covariates
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
      } else {
        # every location has common rho
        rho_i <- rep(rho_sam[i,], nrow(Xpred))
      }
      
      # save the rhos for output
      rho_sample <- rbind(rho_sample,rho_i)
      
      # latent f (linear predictor a + xb + zl)
      # the last zl is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)
      f_sam <- rbind(f_sam,f_i)
      
      # mean for beta distribution using logit-link
      mu_i <- inv_logit(f_i)
      
      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,rho) with common a
      EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
      probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))
      
      # save the expectations and prob. of zeros
      EY_sam <- rbind(EY_sam,EYi)
      probzero_sam <- rbind(probzero_sam, probzero_i)
      
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
      # save the sample
      y_sam <- rbind(y_sam, y_i)
    }
    
    # save predictions for species j to list
    res_list[[sp_name_list[j]]] <- list(f_sam = f_sam,
                                        y_sam = y_sam,
                                        EY_sam = EY_sam,
                                        probzero_sam = probzero_sam,
                                        rho_sample = rho_sample)
  }
  # finally, return list that includes predictions for all J species
  return(res_list)
}

### TEST: does it work?
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
  
  ### sample latent factors for prediction locations
  Zpred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) # (M x n_factor) matrix where M is number of prediction locations
  Z_pi_pred <- matrix(rnorm(nrow(X.pred)*n_factors,0,1),ncol=n_factors) 
  
  # initialize the list including predictions for multiple species
  res_list <- list()
  
  ### loop over species
  n_species <- length(sp_name_list)
  
  for(j in 1:n_species) {
    # load the posterior samples
    alpha_sam <- as.matrix(stan_fit, pars = paste0("alpha[",j,"]"))
    alpha_pi_sam <- as.matrix(stan_fit, pars = paste0("alpha_pi[",j,"]"))
    beta_sam <- as.matrix(stan_fit, pars = c(paste0("beta_1[",1:ncol(X.pred),",",j,"]"),paste0("beta_2[",1:ncol(X.pred),",",j,"]")))
    beta_pi_sam <- as.matrix(stan_fit, pars = c(paste0("beta_pi_1[",1:ncol(X.pred),",",j,"]"),paste0("beta_pi_2[",1:ncol(X.pred),",",j,"]")))
    lambda_sam <- as.matrix(stan_fit, pars = paste0("Lambda[",1:n_factors,",",j,"]"))
    lambda_pi_sam <- as.matrix(stan_fit, pars = paste0("Lambda_pi[",1:n_factors,",",j,"]"))
    
    # prepare rho
    if (rho_modeled) {
      # coefficient related to rho(x)
      alpha_rho_sam <- as.matrix(stan_fit, pars = paste0("alpha_rho[",j,"]"))
      beta_rho_sam <- as.matrix(stan_fit, pars = c(paste0("beta_rho_1[",1:ncol(X.pred),",",j,"]"),paste0("beta_rho_2[",1:ncol(X.pred),",",j,"]")))
    } else {
      # posterior sample of common rho
      rho_sam <- as.matrix(stan_fit, pars = paste0("rho[",j,"]")) #j:th species
    }
    
    # initialize matrices for output
    f_sam <- c()
    y_sam <- c()
    EY_sam <- c()
    probzero_sam <- c()
    prob_suit_sam <- c()
    EY_if_suitable_sam <- c()
    rho_sample <- c()
    
    # loop over posterior samples (take every thinning:th sample)
    for (i in seq(1,nrow(beta_sam),thinning)) {
      # corresponding coefficients
      alpha_i <- alpha_sam[i,]
      alpha_pi_i <- alpha_pi_sam[i,]
      beta_i <- beta_sam[i,]
      beta_pi_i <- beta_pi_sam[i,]
      lambda_i <- lambda_sam[i,]
      lambda_pi_i <- lambda_pi_sam[i,]
      
      # calculate rho
      if (rho_modeled){
        rho_i <- min_rho + (C-min_rho)*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
      } else {
        rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
      }
      
      # save rho for output
      rho_sample <- rbind(rho_sample, rho_i)
      
      # latent f (linear predictor) for mean of beta
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_i <- as.vector(alpha_i + Xpred %*% beta_i + Zpred %*% lambda_i)
      f_sam <- rbind(f_sam,f_i)
      
      # calculate the mean of beta using logit-link
      mu_i <- inv_logit(f_i)
      
      # latent f (linear predictor) for probability of suitability
      # the last ZL is the latent factor * loading that creates interspecific correlations into the model
      f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i + Z_pi_pred %*% lambda_pi_i)
      
      # calculate probability of suitaility using logit-link
      pi_i <- inv_logit(f_pi_i)
      prob_suit_sam <- rbind(prob_suit_sam,pi_i)
      
      ### calculate expectations and probability of zeroes
      ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
      EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
      probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
      
      # save for output
      EY_sam <- rbind(EY_sam,EYi)
      EY_if_suitable_sam <- rbind(EY_if_suitable_sam, EYi_if_suitable)
      probzero_sam <- rbind(probzero_sam, probzero_i)
      
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
      y_sam <- rbind(y_sam, y_i)
    }
    
    # save predictions for species j to list
    res_list[[sp_names[j]]] <- list(f_sam = f_sam,
                                    y_sam = y_sam,
                                    EY_sam = EY_sam,
                                    probzero_sam = probzero_sam,
                                    prob_suit_sam = prob_suit_sam,
                                    EY_if_suitable_sam = EY_if_suitable_sam,
                                    rho_sample = rho_sample)
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

prepare_P_matrix <- function(pred.locations, grid_center.locations) {
  ### prepares a sparse (n_locations X n_spatial_cells) matrix that indicates the spatial random effect block each point belongs to
  ### each rows sums to 1, 1 indicating the random effect block
  ### Given P matrix and vector z of n_spatial_cells spatial grid cells, Pz gives a vector of spatial random effects for each prediction location!
  # pred.locations: includes prediction locations as Nx2-matrix (coordinates in TM35FIN system)
  # grid_center.locations: includes the spatial grid centers as matrix (coordinates in TM35FIN system)
  ### RETURNS matrix P
  
  # turn the matrixes into terra vectors
  # this requires library(terra)
  grid_center.vect <- vect(grid_center.locations, geom = c("x","y"), crs = "EPSG:3067")
  pred_locs.vect <- vect(pred.locations, geom = c("x","y"), crs = "EPSG:3067")
  
  # find nearest spatial grid cell for each prediction location
  nearest_grid_cells.df <- as.data.frame(nearest(pred_locs.vect, grid_center.vect))
  # vector telling the ID of the closest grid cell
  nearest_id <- nearest_grid_cells.df$to_id 
  
  # initialize P as zeros
  P <- matrix(0, nrow = nrow(pred.locations), ncol = nrow(grid_center.locations))
  # loop over prediction locations, add 1 to the column corresponding to the closest spatial grid cell
  for (i in 1:nrow(pred.locations)) {
    P[i,nearest_id[i]] <- 1
  }
  
  # name rows and columns & return
  colnames(P) <- 1:nrow(grid_center.locations)
  rownames(P) <- 1:nrow(pred.locations)
  return(P)
}


predict_spatial_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, sp_name_list, n_factors, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000) {
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
  

  ### sample latent factors for prediction locations
  for (k in 1:n_factors) {
    break
  }
  
  
  # load in posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  
  # prepare rho
  if (rho_modeled) {
    # rho(x) modeled with covariates
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # every location has common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }
  
  # spatial random effects
  phi_sam <- as.matrix(stan_fit, pars = c("phi"))
  # covariance function parameters
  l_sam <- as.matrix(stan_fit, pars = c("l"))
  s2_sam <- as.matrix(stan_fit, pars = c("s2_cf"))
  
  
  # initialize matrices for output
  f_sam <- c()
  phi_pred_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    beta_i <- beta_sam[i,]
    
    # calculate rho
    if (rho_modeled){
      rho_i <- C*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
    }
    
    # save for output
    rho_sample <- rbind(rho_sample, rho_i)
    
    # latent f (linear predictor)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    
    # sample spatial random effect
    phi_obs_i <- phi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_sam[i,]
    s2_i <- s2_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### sample spatial random effects over the prediction locations
    # add jitter on diagonal for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial random effects
    phi_pred_sam <- rbind(phi_pred_sam, as.vector(P %*% phi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_i <- f_i + as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
    f_sam <- rbind(f_sam,f_i)
    
    # mean of the beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,rho) with common a
    EYi <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a))
    probzero_i <- mapply(calculate_probzero_LC_beta, mu=mu_i, rho=rho_i, MoreArgs = list(a=a))
    
    # save the expectations and prob. of zeros
    EY_sam <- rbind(EY_sam,EYi)
    probzero_sam <- rbind(probzero_sam, probzero_i)
    
    ### also predict Ys
    # latent beta variables
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    # scale & censor
    y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              phi_mu_pred_sam = phi_pred_sam,
              rho_sample = rho_sample))
}

predict_spatial_ZI_beta_regression_JSDM <- function(stan_fit, X.pred, X.orig, pred.locations, S.pred, S.obs, thinning = 10, a = 1, rho_modeled = FALSE, C = 1000) {
  ### function to make predictions with zero-inflated spatial left-censored beta regression
  # X.pred: prediction matrix (mxp)
  # X.orig: non-scaled data matrix (nxp) to learn about the scaling parameters
  # pred.locations: locations for each prediction cell (in meters)
  # S.pred: locations of coarse spatial grid cells (in kilometers)
  # S.obs: locations of observed coarse spatial grid cells (in kilometers, WATCH THAT THIS IS THE SAME AS WHEN FITTING A MODEL)
  # thinning: use every thinning:th posterior sample in predictions
  # a: left-censoring constant
  # rho_modeled: TRUE: rho modeled with covariates, FALSE: common rho
  # C: if rho modeled, what is it's upper limit? rho = C*inv_logit(a+XB)
  # RETURNS list of rep x m matrix of samples from posterior predictive for different quantities
  ### 1) latent linear predictors f 
  ### 2) predicted Ys
  ### 3) expected Ys
  ### 4) probabilies of zero
  ### 5) probability of suitability
  ### 6) predicted spatial random effects (for mean of beta distribution)
  ### 7) predicted spatial random effects (for probability of suitability)
  ### 8) rhos
  
  # load in posterior samples
  alpha_sam <- as.matrix(stan_fit, pars = c("alpha"))
  alpha_pi_sam <- as.matrix(stan_fit, pars = c("alpha_pi"))
  beta_sam <- as.matrix(stan_fit, pars = c("beta_1","beta_2"))
  beta_pi_sam <- as.matrix(stan_fit, pars = c("beta_pi_1","beta_pi_2"))
  
  # prepare rho
  if (rho_modeled) {
    # rho(x) modeled with covariates
    alpha_rho_sam <- as.matrix(stan_fit, pars = c("alpha_rho"))
    beta_rho_sam <- as.matrix(stan_fit, pars = c("beta_rho_1","beta_rho_2"))
  } else {
    # every location has common rho
    rho_sam <- as.matrix(stan_fit, pars = c("rho"))
  }
  
  # spatial random effects (mean of beta & prob. of suitability)
  phi_sam <- as.matrix(stan_fit, pars = c("phi_mu"))
  phi_pi_sam <- as.matrix(stan_fit, pars = c("phi_pi"))
  
  # covariance function parameters (for mean of beta)
  l_sam <- as.matrix(stan_fit, pars = c("l_mu"))
  s2_sam <- as.matrix(stan_fit, pars = c("s2_mu"))
  
  # covariance function parameters (for prob. of suitability)
  l_pi_sam <- as.matrix(stan_fit, pars = c("l_pi"))
  s2_pi_sam <- as.matrix(stan_fit, pars = c("s2_pi"))
  
  # prepare prediction matrix (scale, add second order terms)
  X.pred.scaled <- scale_covariates(X.orig,X.pred)
  Xpred <- add_second_order_terms(X.pred.scaled, colnames(X.pred))
  Xpred <- as.matrix(Xpred) # from df to matrix for calculation
  
  # prepare P matrix (mxS) that tells in which spatial random effect cell each prediction point belongs to
  P <- prepare_P_matrix(pred.locations,S.pred*1000) # turn grid locations to meters
  
  # initialize matrices for output
  f_sam <- c()
  y_sam <- c()
  EY_sam <- c()
  probzero_sam <- c()
  prob_suit_sam <- c()
  EY_if_suitable_sam <- c()
  phi_mu_pred_sam <- c()
  phi_pi_pred_sam <- c()
  rho_sample <- c()
  
  # loop over posterior samples (take every thinning:th sample)
  for (i in seq(1,nrow(beta_sam),thinning)) {
    # corresponding coefficients
    alpha_i <- alpha_sam[i,]
    alpha_pi_i <- alpha_pi_sam[i,]
    beta_i <- beta_sam[i,]
    beta_pi_i <- beta_pi_sam[i,]
    
    # calculate rho
    if (rho_modeled){
      rho_i <- C*inv_logit(as.vector(alpha_rho_sam[i,]  + Xpred %*% beta_rho_sam[i,]))
    } else {
      rho_i <- rep(rho_sam[i,], nrow(Xpred)) # every location has common rho
    }
    
    # save for output
    rho_sample <- rbind(rho_sample, rho_i)
    
    ### 1) PREDICT SPATIAL RANDOM EFFECTS FOR MEAN OF BETA DISTRIBUTION
    # latent f (linear predictor)
    f_i <- as.vector(alpha_i + Xpred %*% beta_i)
    
    phi_obs_i <- phi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_sam[i,]
    s2_i <- s2_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### add jitter for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial effects
    phi_mu_pred_sam <- rbind(phi_mu_pred_sam, as.vector(P %*% phi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_i <- f_i + as.vector(P %*% phi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
    f_sam <- rbind(f_sam,f_i)
    
    # mean of the beta distribution using logit-link
    mu_i <- inv_logit(f_i)
    
    ### 2) PREDICT SPATIAL RANDOM EFFECTS FOR PROBABILITY OF ZERO
    # latent f (linear predictor)
    f_pi_i <- as.vector(alpha_pi_i + Xpred %*% beta_pi_i)
    
    # sample spatial random effects
    phi_obs_i <- phi_pi_sam[i,] #these are spatial random effects that are sampled (included in the training data), correspond to locations S.obs
    l_i <- l_pi_sam[i,]
    s2_i <- s2_pi_sam[i,]
    
    # calculate covariance block-matrices
    K_pred_obs <- exp_covariance(S.pred,S.obs,s2_i,l_i)
    K_pred <- exp_covariance(S.pred,S.pred,s2_i,l_i)
    K_obs <- exp_covariance(S.obs,S.obs,s2_i,l_i)
    
    # mean and covariance for predicting phi in new locations given set of parameters and observed phi
    # justification for formulas can be read from Pietilä thesis (2025)
    phi_pred_m <- K_pred_obs %*% solve(K_obs,phi_obs_i)
    phi_pred_Cov <- K_pred - K_pred_obs %*% solve(K_obs) %*% t(K_pred_obs)
    
    ### add jitter for computational stability
    phi_pred_Cov <- phi_pred_Cov + 1e-08*diag(length(phi_pred_m))
    
    # use Cholesky for sampling (if Sigma = LL^t and I ~ MVN(0,1), then u + LI ~ MVN(u,Sigma))
    L <- t(chol(phi_pred_Cov))
    phi_pi_pred_i <- phi_pred_m + L %*% rnorm(length(phi_pred_m))
    
    # save spatial effects
    phi_pi_pred_sam <- rbind(phi_pi_pred_sam, as.vector(P %*% phi_pi_pred_i))
    
    # add the spatial random effect to linear predictor
    f_pi_i <- f_pi_i + as.vector(P %*% phi_pi_pred_i) # P matrix picks the correct random effects (from the grid cell that the prediction point belongs to)
    
    ### probability of suitability using logit-link
    pi_i <- inv_logit(f_pi_i)
    prob_suit_sam <- rbind(prob_suit_sam, pi_i)
    
    ### calculate expectations and probability of zeroes
    ### NOTE: mapply iterates over pairs of (mu,pi,rho) with common a
    EYi <- mapply(integrate_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    EYi_if_suitable <- mapply(integrate_LC_beta,mu=mu_i,rho=rho_i, MoreArgs = list(a=a)) #NOTE: this does not care about prob. of suitability
    probzero_i <- mapply(calculate_probzero_ZI_LC_beta,mu=mu_i,pi=pi_i,rho=rho_i, MoreArgs = list(a=a))
    
    # save for output
    EY_sam <- rbind(EY_sam,EYi)
    EY_if_suitable_sam <- rbind(EY_if_suitable_sam, EYi_if_suitable)
    probzero_sam <- rbind(probzero_sam, probzero_i)
    
    ### also predict Ys
    # sample suitabilities of locations
    Z_i <- rbinom(nrow(Xpred),1,pi_i) #0/1 vector to tell whether the location is suitable (1 => Y from beta) or unsuitable (0 => Y=0)
    # sample latent beta variables
    V_i <- rbeta(nrow(Xpred),mu_i*rho_i,(1-mu_i)*rho_i)
    # scale & censor
    y_i <- sapply(V_i, function(v) (max(0,(a+1)*v - a)))
    # if the location was unsuitable, automatically y = 0
    y_i <- Z_i*y_i 
    y_sam <- rbind(y_sam, y_i)
  }
  
  # return everything
  return(list(f_sam = f_sam,
              y_sam = y_sam,
              EY_sam = EY_sam,
              probzero_sam = probzero_sam,
              prob_suit_sam = prob_suit_sam,
              EY_if_suitable_sam = EY_if_suitable_sam,
              phi_mu_pred_sam = phi_mu_pred_sam,
              phi_pi_pred_sam = phi_pi_pred_sam,
              rho_sample = rho_sample))
}